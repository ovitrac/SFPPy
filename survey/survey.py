"""
survey/survey.py — Main Survey Class for Exposure Estimation
============================================================

Production-grade Survey class for survey-scale migration estimation.

Features:
- Parallel master curve computation with persistent cache
- Dynamic PDF generation from deterministic priors
- Add/remove substances with cache reuse
- Resume support for interrupted computations
- Multilayer support with reference layer selection

@project: SFPPy — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

import uuid
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple, Union
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict

import numpy as np

from survey.models import (
    LayerSpec,
    PackagingSpec,
    SubstanceSpec,
    PriorSpec,
    SurveyConfig,
    SubstanceModel,
    select_reference_layer,
)
from survey.priors import discretize_prior, compute_weight_matrix
from survey.fingerprints import (
    fingerprint_physics,
    fingerprint_probability,
    fingerprint_survey_config,
)
from survey.cache import (
    MasterCurveCache,
    MasterCurveKey,
    PdfCache,
    SurveyState,
    SurveyStateManager,
)
from survey.workers import worker_compute_curve
from survey.reporting import ProgressReporter, format_survey_preview, format_pdf_summary
from survey.io import load_scenario, load_substances_from_scenario, save_results, save_manifest


class Survey:
    """
    Survey-scale migration estimation with master curve caching.

    The Survey class provides:
    - Parallel computation of substance-specific master curves
    - Persistent caching to avoid redundant computation
    - Dynamic PDF computation from priors
    - Add/remove substances with minimal recomputation
    - Resume capability for interrupted runs
    - Multilayer support with reference layer selection

    Example
    -------
    >>> from survey import Survey
    >>> survey = Survey.from_scenario("scenario.yml")
    >>> survey.compute()
    >>> print(survey.summary())
    >>> survey.save("output.npz")
    """

    def __init__(
        self,
        config: SurveyConfig,
        substances: List[SubstanceSpec],
        substance_model: Optional[SubstanceModel] = None,
    ):
        """
        Initialize Survey.

        Parameters
        ----------
        config : SurveyConfig
            Survey configuration.
        substances : List[SubstanceSpec]
            List of substances to include.
        substance_model : SubstanceModel, optional
            Model for parameter inference. If None, uses default.
        """
        self.config = config
        self._substances: List[SubstanceSpec] = []
        self._substance_model = substance_model or SubstanceModel(
            food_simulant=config.packaging.food_simulant
        )

        # Caches
        self._curve_cache = MasterCurveCache(config.cache_dir)
        self._pdf_cache = PdfCache(config.cache_dir)
        self._state_manager = SurveyStateManager(config.cache_dir)

        # State
        self._survey_id = str(uuid.uuid4())[:8]
        self._i_ref = 0
        self._ref_layer_details: List[Dict] = []
        self._curves: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
        self._results: Optional[Dict[str, np.ndarray]] = None

        # Select reference layer (for multilayer)
        self._select_reference_layer()

        # Add substances (triggers inference)
        for s in substances:
            self.add_substance(s)

    @classmethod
    def from_scenario(cls, path: Union[str, Path]) -> "Survey":
        """
        Create Survey from scenario YAML file.

        Parameters
        ----------
        path : str or Path
            Path to scenario YAML.

        Returns
        -------
        Survey
            Initialized Survey instance.
        """
        config = load_scenario(path)
        substances = load_substances_from_scenario(path)
        return cls(config=config, substances=substances)

    def _select_reference_layer(self) -> None:
        """Select reference layer for multilayer packaging."""
        layers = self.config.packaging.layers
        if len(layers) == 1:
            self._i_ref = 0
            self._ref_layer_details = [{
                'index': 0,
                'polymer': layers[0].polymer,
                'thickness_m': layers[0].thickness_m,
            }]
        else:
            self._i_ref, self._ref_layer_details = select_reference_layer(
                layers, self._substance_model
            )

    @property
    def ref_layer(self) -> LayerSpec:
        """Reference layer for Fo definition."""
        return self.config.packaging.layers[self._i_ref]

    @property
    def substances(self) -> List[SubstanceSpec]:
        """List of substances with inferred parameters."""
        return list(self._substances)

    @property
    def n_substances(self) -> int:
        """Number of substances."""
        return len(self._substances)

    def add_substance(self, substance: SubstanceSpec) -> None:
        """
        Add substance to survey.

        Triggers parameter inference and uniqueness check.

        Parameters
        ----------
        substance : SubstanceSpec
            Substance to add.

        Raises
        ------
        ValueError
            If duplicate substance (by canonical ID).
        """
        # Check uniqueness
        new_id = substance.canonical_id()
        for existing in self._substances:
            if existing.canonical_id() == new_id:
                raise ValueError(f"Duplicate substance: {new_id}")

        # Infer parameters using reference layer
        inferred = self._substance_model.infer_all(substance, self.ref_layer)
        self._substances.append(inferred)

        # Invalidate results (require recomputation)
        self._results = None

    def remove_substance(self, substance_id: str) -> bool:
        """
        Remove substance by ID.

        Parameters
        ----------
        substance_id : str
            Substance ID or canonical ID to remove.

        Returns
        -------
        bool
            True if removed, False if not found.
        """
        for i, s in enumerate(self._substances):
            if s.id == substance_id or s.canonical_id() == substance_id:
                del self._substances[i]
                self._results = None
                return True
        return False

    def _build_curve_tasks(self) -> List[Dict[str, Any]]:
        """Build task payloads for master curve computation."""
        tasks = []
        layer = self.ref_layer
        pkg = self.config.packaging

        # Discretize time prior to get max time
        time_vals, _ = discretize_prior(self.config.time_prior)
        t_max = float(np.max(time_vals))
        l2 = layer.thickness_m ** 2

        for substance in self._substances:
            D = substance.D or 1e-14
            k = substance.k or 1.0
            k0 = substance.k0 or 1.0

            Fo_max = (D * t_max / l2) * self.config.fo_max_factor

            key = MasterCurveKey(
                polymer=layer.polymer,
                mass_g_mol=substance.mass_g_mol,
                D=D,
                k=k,
                k0=k0,
                lP_m=layer.thickness_m,
                Fo_max=float(Fo_max),
                h=pkg.h_m_s,
                surface_area=pkg.surface_area_m2,
                food_volume=pkg.food_volume_m3,
                contact_temperature_degC=pkg.contact_temperature_degC,
                CF0=pkg.cf0,
                n_fo=self.config.n_fo,
            )

            tasks.append({
                "cache_dir": self.config.cache_dir,
                "key_dict": key.to_dict(),
                "substance_id": substance.id,
                "fo_min_floor": self.config.fo_min_floor,
            })

        return tasks

    def build_master_curves(
        self,
        parallel: bool = True,
        max_workers: Optional[int] = None,
    ) -> Dict[str, int]:
        """
        Build or load master curves for all substances.

        Parameters
        ----------
        parallel : bool
            Use parallel workers (default: True).
        max_workers : int, optional
            Maximum parallel workers. If None, uses CPU count.

        Returns
        -------
        dict
            Statistics: {'hits': int, 'misses': int}.
        """
        tasks = self._build_curve_tasks()
        if not tasks:
            return {'hits': 0, 'misses': 0}

        # Pre-check cache for progress display
        self._curve_cache.reset_stats()
        prog = ProgressReporter(total_steps=len(tasks), desc="Master Curves")

        for task in tasks:
            key = MasterCurveKey(**task["key_dict"])
            if self._curve_cache.exists(key):
                prog.hit()
            else:
                prog.miss()

        # Compute curves
        if parallel and len(tasks) > 1:
            with ProcessPoolExecutor(max_workers=max_workers) as exe:
                futures = {
                    exe.submit(worker_compute_curve, task): task["substance_id"]
                    for task in tasks
                }
                for fut in as_completed(futures):
                    result = fut.result()
                    substance_id = futures[fut]
                    self._curves[substance_id] = (result["fo"], result["cf"])
                    prog.update()
        else:
            for task in tasks:
                result = worker_compute_curve(task)
                self._curves[task["substance_id"]] = (result["fo"], result["cf"])
                prog.update()

        prog.done()

        return self._curve_cache.stats

    def _family_mass_ceiling(self, per_substance_ceiling: float) -> float:
        """Upper bound on the combined family ``CF/CP0`` given that every
        substance's own plateau is ``<= per_substance_ceiling`` (= V_focal/V_F).

        Mirrors the probabilistic family-combination model exactly:
        - **exchangeable** substances are alternatives — their combined
          contribution is the occurrence-weighted **average** $\\sum_j p_j g_j$
          ($\\sum_j p_j = 1$), a convex combination of the per-substance curves,
          so it is bounded by the largest member and ``<= per_substance_ceiling``
          (one alternative present at the family concentration; no 1/E — F8,
          2026-06-12, consistent with the kernel's exchangeable term);
        - **non-exchangeable** substances are co-present and contribute their
          normalised-weight sum, ``<= per_substance_ceiling × Σ (w_j/max w)``.

        Used by the family-level mass-conservation check so that a legitimate
        multi-substance family (whose total migration is the weighted sum over
        co-present substances) is not falsely flagged against a single-substance
        bound, while a genuine over-prediction still trips it. The per-substance
        guard in the solver (``survey.workers.assert_bilayer_mass_conservation``)
        remains the independent first check.
        """
        subs = self._substances
        exch = [s for s in subs if getattr(s, 'exchangeable', True)]
        non = [s for s in subs if not getattr(s, 'exchangeable', True)]
        total = 0.0
        if exch:
            # weighted average of alternatives (Σ p_j = 1) ≤ max member ≤ ceiling
            total += per_substance_ceiling
        if non:
            w = np.array([getattr(s, 'weight', 1.0) for s in non], dtype=float)
            w_max = w.max()
            w_norm = (w / w_max) if w_max > 0 else np.ones(len(non))
            total += per_substance_ceiling * float(w_norm.sum())
        return total

    def _compute_cf_tensor(self) -> Dict[str, np.ndarray]:
        """
        Compute CF tensor from master curves and priors.

        Handles exchangeability of substances:
        - Exchangeable: Alternatives (one OR another is present).
          Concentration is fractionated among exchangeable substances.
          Represents uncertainty about which specific substance is used.
        - Non-exchangeable: Always present (AND).
          Uses full family concentration with normalized weights.
          Monomers are non-exchangeable by definition.

        Mathematical formulation:
        - Exchangeable: CF_exch = g_exch(t) × C0/E where E = number of exchangeable
        - Non-exchangeable: CF_non = Σ (w_j/max(w)) × g_j(t) × C0
        - Combined: CF_family = CF_exch + CF_non
        """
        # Discretize priors
        time_vals, time_w = discretize_prior(self.config.time_prior)
        conc_vals, conc_w = discretize_prior(self.config.conc_prior)

        layer = self.ref_layer
        l2 = layer.thickness_m ** 2

        # Build D array for all substances
        D_vals = np.array([s.D or 1e-14 for s in self._substances], dtype=float)
        k_vals = np.array([s.k or 1.0 for s in self._substances], dtype=float)
        k0_vals = np.array([s.k0 or 1.0 for s in self._substances], dtype=float)

        n_sub = len(self._substances)
        n_t = len(time_vals)

        # Fo matrix: (N_T, n_sub)
        Fo = np.outer(time_vals, D_vals / l2)

        # Interpolate master curves at Fo values
        g_tj = np.empty((n_t, n_sub), dtype=float)
        for j, substance in enumerate(self._substances):
            fo_grid, cf_over_cp0 = self._curves[substance.id]
            g_tj[:, j] = np.interp(
                np.clip(Fo[:, j], 0.0, None),
                fo_grid,
                cf_over_cp0,
            )

        # Separate substances by exchangeability
        # BACKWARD COMPATIBILITY: If 'exchangeable' attribute is missing, default to True
        # This ensures older scenarios without the field work as before (all substances
        # treated as exchangeable alternatives with concentration fractionation)
        exch_indices = [i for i, s in enumerate(self._substances)
                        if getattr(s, 'exchangeable', True)]
        non_exch_indices = [i for i, s in enumerate(self._substances)
                           if not getattr(s, 'exchangeable', True)]

        n_exch = len(exch_indices)
        n_non_exch = len(non_exch_indices)

        # Initialize combined migration curve
        g_combined = np.zeros(n_t, dtype=float)

        # Track contribution fractions for diagnostics
        conc_fraction_exch = 0.0
        conc_fraction_non_exch = 0.0

        # =========================================================
        # Exchangeable substances: mixture model with fractionation
        # =========================================================
        if n_exch > 0:
            # Get occurrence weights for exchangeable substances
            w_exch = np.array([
                getattr(self._substances[i], 'weight', 1.0)
                for i in exch_indices
            ], dtype=float)

            # Normalize to probabilities (sum to 1)
            w_sum = w_exch.sum()
            if w_sum > 0:
                p_exch = w_exch / w_sum
            else:
                p_exch = np.ones(n_exch) / n_exch

            # Weighted average migration curve for exchangeable
            g_exch_curves = g_tj[:, exch_indices]  # Shape: (n_t, n_exch)
            g_exch = np.dot(g_exch_curves, p_exch)  # Shape: (n_t,)

            # Exchangeable contribution: weighted average of the alternatives.
            # CF_exch = C0 * Σ_j p_j g_j  (one alternative present at the family
            # concentration C0). The former 1/E factor was a SECOND dilution on
            # top of the already-normalised weighted mean (Σ p_j = 1) and is
            # removed (decision 2026-06-12; see REVIEW_20260612*, finding F8).
            conc_fraction_exch = 1.0
            g_combined += g_exch * conc_fraction_exch

        # =========================================================
        # Non-exchangeable substances: full concentration, normalized weights
        # Each substance contributes independently: CF_j = g_j(t) × C0 × w_norm_j
        # Example: styrene (occ=2) and EB (occ=1) at 300 ppm
        #   w_norm = [2/2, 1/2] = [1.0, 0.5]
        #   styrene: 300×1.0 = 300 ppm, EB: 300×0.5 = 150 ppm
        # =========================================================
        w_norm = np.array([])  # Initialize for scope
        if n_non_exch > 0:
            # Get occurrence weights for non-exchangeable substances
            w_non = np.array([
                getattr(self._substances[i], 'weight', 1.0)
                for i in non_exch_indices
            ], dtype=float)

            # Normalize so maximum weight is 1 (not sum to 1)
            # This means the highest-occurrence non-exchangeable gets full concentration
            w_max = w_non.max()
            if w_max > 0:
                w_norm = w_non / w_max
            else:
                w_norm = np.ones(n_non_exch)

            # Each non-exchangeable substance contributes independently
            # Sum of weighted migration curves: Σ w_norm_j × g_j(t)
            # This is NOT averaged - each substance adds its own contribution
            g_non_curves = g_tj[:, non_exch_indices]  # Shape: (n_t, n_non_exch)
            g_non = np.dot(g_non_curves, w_norm)  # Shape: (n_t,) - weighted sum, not average

            # Non-exchangeable contribution: each substance contributes w_norm_j × C0
            # Total contribution = Σ g_j × C0 × w_norm_j = C0 × Σ (g_j × w_norm_j) = C0 × g_non
            conc_fraction_non_exch = 1.0  # Full concentration, weighting is in g_non
            g_combined += g_non * conc_fraction_non_exch

        # Store sum of normalized weights for diagnostics
        # For non-exch with styrene(2)+EB(1): w_norm_sum = 1.0 + 0.5 = 1.5
        w_norm_sum = w_norm.sum() if n_non_exch > 0 else 0.0

        # If no substances at all, fall back to zero migration
        if n_sub == 0:
            g_combined = np.zeros(n_t)
            conc_fraction_exch = 0.0
            conc_fraction_non_exch = 0.0
            w_norm_sum = 0.0

        # CF_family(t_i, CP0_j) = g_combined(t_i) * CP0_j
        cf_family = np.outer(g_combined, conc_vals)

        # Family-level mass-conservation check (proper probabilistic weights).
        # Per-substance conservation (g_j <= V_focal/V_F) is the solver's job;
        # here we bound the COMBINED family curve by the SAME exchangeable/
        # non-exchangeable weighting applied to that per-substance ceiling, so a
        # legitimate multi-substance family is not flagged against a single-
        # substance bound. Monolayer: focal layer == the only layer.
        import warnings
        from survey.workers import MASS_CONSERVATION_TOL
        V_F = self.config.packaging.food_volume_m3
        if V_F > 0:
            per_sub_ceiling = layer.thickness_m * \
                self.config.packaging.surface_area_m2 / V_F
            family_ceiling = self._family_mass_ceiling(per_sub_ceiling)
            g_max = float(g_combined.max())
            if g_max > family_ceiling * (1.0 + MASS_CONSERVATION_TOL):
                warnings.warn(
                    f"family CF/CP0 max={g_max:.4g} exceeds weighted family mass "
                    f"bound={family_ceiling:.4g} (per-substance V_focal/V_F="
                    f"{per_sub_ceiling:.4g}); check substance combination weights",
                    RuntimeWarning,
                )

        # Build substance weights array for backward compatibility
        substance_weights = np.array([
            getattr(s, 'weight', 1.0) for s in self._substances
        ], dtype=float)
        w_sum = substance_weights.sum()
        if w_sum > 0:
            substance_weights /= w_sum
        else:
            substance_weights = np.ones(n_sub) / max(n_sub, 1)

        return {
            'CF_tensor': cf_family,
            'CF_samples': cf_family.ravel(),
            'weights': compute_weight_matrix(time_w, conc_w).ravel(),
            'Fo': Fo,
            'D_vals': D_vals,
            'k_vals': k_vals,
            'k0_vals': k0_vals,
            'conc_vals': conc_vals,
            'conc_weights': conc_w,
            'time_vals': time_vals,
            'time_weights': time_w,
            'substance_weights': substance_weights,  # Occurrence-based weights
            # Exchangeability diagnostics
            'n_exchangeable': n_exch,
            'n_non_exchangeable': n_non_exch,
            'conc_fraction_exchangeable': conc_fraction_exch,
            'conc_fraction_non_exchangeable': conc_fraction_non_exch,
            'non_exch_weight_sum': w_norm_sum,  # Sum of normalized weights for non-exchangeable
        }

    def _per_substance_factors(self) -> np.ndarray:
        """Per-substance contribution factor: exchangeable → p_i = w_i/Σ_exch w
        (expected presence); non-exchangeable → w_i/max_non w (co-presence).
        Such that Σ_i factor_i · g_i == the combined family curve/surface."""
        subs = self._substances
        n_sub = len(subs)
        factor = np.zeros(n_sub, dtype=float)
        exch_idx = [i for i, s in enumerate(subs) if getattr(s, 'exchangeable', True)]
        non_idx = [i for i, s in enumerate(subs) if not getattr(s, 'exchangeable', True)]
        if exch_idx:
            w = np.array([getattr(subs[i], 'weight', 1.0) for i in exch_idx], dtype=float)
            ws = w.sum()
            p = (w / ws) if ws > 0 else np.ones(len(exch_idx)) / len(exch_idx)
            for k, i in enumerate(exch_idx):
                factor[i] = p[k]
        if non_idx:
            w = np.array([getattr(subs[i], 'weight', 1.0) for i in non_idx], dtype=float)
            wm = w.max()
            wn = (w / wm) if wm > 0 else np.ones(len(non_idx))
            for k, i in enumerate(non_idx):
                factor[i] = wn[k]
        return factor

    def per_substance_cf_tensors(self) -> Dict[str, Dict[str, np.ndarray]]:
        """Decompose the family CF into per-substance contributions.

        Concentration is applied AFTER the (C0-independent) master curve, so
        each substance's *expected* contribution is a post-hoc re-weighting of
        its own curve g_i:

          - exchangeable     : CF_i = p_i · outer(g_i, conc_vals),  p_i = w_i/Σ_exch w
          - non-exchangeable : CF_i = (w_i/max_non w) · outer(g_i, conc_vals)

        By construction ``Σ_i CF_i == family CF_tensor`` (drop-1/E model). Used
        by the per-substance aggregation (keyed by CAS) requested by Mai; needs
        no re-solve — only the cached master curves. Base Survey = monolayer /
        single-curve path; BilayerSurvey overrides for its surface combination.

        Returns
        -------
        dict  {cas: {CF_tensor, CF_samples, weights, time_vals, time_weights,
                     conc_vals, conc_weights, exchangeable, weight}}
        """
        time_vals, time_w = discretize_prior(self.config.time_prior)
        conc_vals, conc_w = discretize_prior(self.config.conc_prior)
        layer = self.ref_layer
        l2 = layer.thickness_m ** 2
        subs = self._substances
        n_sub = len(subs)
        D_vals = np.array([s.D or 1e-14 for s in subs], dtype=float)
        Fo = np.outer(time_vals, D_vals / l2)

        g_tj = np.empty((len(time_vals), n_sub), dtype=float)
        for j, s in enumerate(subs):
            fo_grid, cf_over_cp0 = self._curves[s.id]
            g_tj[:, j] = np.interp(np.clip(Fo[:, j], 0.0, None), fo_grid, cf_over_cp0)

        # Per-substance weight factor (same convention as _compute_cf_tensor).
        factor = self._per_substance_factors()

        wmat = compute_weight_matrix(time_w, conc_w).ravel()
        out: Dict[str, Dict[str, np.ndarray]] = {}
        for j, s in enumerate(subs):
            cf_j = np.outer(factor[j] * g_tj[:, j], conc_vals)
            cas = getattr(s, 'cas', None) or getattr(s, 'id', f"idx{j}")
            if cas in out:                       # same CAS twice in a family: accumulate
                out[cas]['CF_tensor'] = out[cas]['CF_tensor'] + cf_j
                out[cas]['CF_samples'] = out[cas]['CF_tensor'].ravel()
                continue
            out[cas] = {
                'CF_tensor': cf_j, 'CF_samples': cf_j.ravel(), 'weights': wmat,
                'time_vals': time_vals, 'time_weights': time_w,
                'conc_vals': conc_vals, 'conc_weights': conc_w,
                'exchangeable': bool(getattr(s, 'exchangeable', True)),
                'weight': float(getattr(s, 'weight', 1.0)),
            }
        return out

    def _compute_pdf(
        self,
        samples: np.ndarray,
        weights: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute weighted PDF and CDF."""
        w = np.asarray(weights, dtype=float)
        s = w.sum()
        if s <= 0:
            raise ValueError("Sum of weights must be positive")
        w /= s

        # Degenerate distributions: when the sample range is zero or denormal
        # (e.g. a non-migrating overpack chain, max CF ~1e-310), np.histogram
        # with integer bins raises "Too many bins for data range. Cannot
        # create N finite-sized bins" (caught on a degenerate near-zero
        # overpack case). Values below 1e-250 are hundreds of decades under
        # LOD/10 — physically zero: emit a single-spike PDF at 0 on [0, 1].
        smp = np.asarray(samples, dtype=float)
        smax = float(np.nanmax(smp)) if smp.size else 0.0
        smin = float(np.nanmin(smp)) if smp.size else 0.0
        degenerate = (not np.isfinite(smax)) or smax < 1e-250 \
            or (smax - smin) <= 0.0 or not np.isfinite(smax - smin)
        if degenerate:
            # physically-zero case → spike in the first bin of [0, 1];
            # equal-but-finite case → spike at the common value.
            if (not np.isfinite(smax)) or smax < 1e-250:
                lo, hi = 0.0, 1.0
            else:
                lo, hi = 0.9 * smax, 1.1 * smax
            edges = np.linspace(lo, hi, self.config.pdf_bins + 1)
            centers = 0.5 * (edges[:-1] + edges[1:])
            hist = np.zeros(self.config.pdf_bins)
            j = self.config.pdf_bins // 2 if lo > 0.0 else 0
            hist[j] = 1.0 / (edges[1] - edges[0])   # all mass in one bin
            cdf = np.concatenate([np.zeros(j), np.ones(self.config.pdf_bins - j)])
            return centers, hist, cdf

        hist, edges = np.histogram(
            samples,
            bins=self.config.pdf_bins,
            weights=w,
            density=True,
        )
        centers = 0.5 * (edges[:-1] + edges[1:])
        cdf = np.cumsum(hist * np.diff(edges))
        cdf = np.minimum(cdf, 1.0)

        return centers, hist, cdf

    def compute(self, parallel: bool = True, max_workers: Optional[int] = None) -> None:
        """
        Compute survey results.

        Builds master curves (if needed) and computes PDF.

        Parameters
        ----------
        parallel : bool
            Use parallel workers for curve computation.
        max_workers : int, optional
            Maximum parallel workers.
        """
        # Build master curves
        self.build_master_curves(parallel=parallel, max_workers=max_workers)

        # Compute CF tensor
        tensor_results = self._compute_cf_tensor()

        # Compute PDF
        centers, pdf, cdf = self._compute_pdf(
            tensor_results['CF_samples'],
            tensor_results['weights'],
        )

        # Store results
        self._results = {
            **tensor_results,
            'pdf_bin_centers': centers,
            'pdf': pdf,
            'cdf': cdf,
        }

    @property
    def results(self) -> Dict[str, np.ndarray]:
        """
        Get computed results.

        Returns
        -------
        dict
            Results dictionary with PDF, CDF, and tensor outputs.

        Raises
        ------
        RuntimeError
            If compute() has not been called.
        """
        if self._results is None:
            raise RuntimeError("Call compute() first")
        return self._results

    def quantile(self, q: float) -> float:
        """
        Get quantile value from CDF.

        Parameters
        ----------
        q : float
            Quantile (0-1).

        Returns
        -------
        float
            Value at quantile.
        """
        if self._results is None:
            raise RuntimeError("Call compute() first")

        cdf = self._results['cdf']
        centers = self._results['pdf_bin_centers']
        idx = np.searchsorted(cdf, q)
        if idx >= len(centers):
            return centers[-1]
        return centers[idx]

    def summary(self) -> str:
        """
        Get summary of computed results.

        Returns
        -------
        str
            Formatted summary string.
        """
        if self._results is None:
            return "No results computed. Call compute() first."

        q50 = self.quantile(0.50)
        q95 = self.quantile(0.95)
        q99 = self.quantile(0.99)
        mean = np.average(
            self._results['CF_samples'],
            weights=self._results['weights'],
        )

        return format_pdf_summary(
            pdf_bins=self.config.pdf_bins,
            q50=q50,
            q95=q95,
            q99=q99,
            mean=mean,
        )

    def preview(self) -> str:
        """
        Get survey preview showing configuration.

        Returns
        -------
        str
            Formatted preview string.
        """
        time_vals, _ = discretize_prior(self.config.time_prior)
        t_max = float(np.max(time_vals))
        l2 = self.ref_layer.thickness_m ** 2

        fo_min = 0.0
        fo_max = 0.0
        for s in self._substances:
            D = s.D or 1e-14
            fo_max = max(fo_max, (D * t_max / l2) * self.config.fo_max_factor)

        return format_survey_preview(
            name=self.config.name,
            n_substances=self.n_substances,
            n_layers=self.config.packaging.n_layers,
            i_ref=self._i_ref,
            polymer=self.ref_layer.polymer,
            thickness_m=self.config.packaging.total_thickness_m,
            fo_range=(fo_min, fo_max),
            cache_dir=self.config.cache_dir,
            cache_entries=len(self._curve_cache.list_entries()),
        )

    def __repr__(self) -> str:
        return (
            f"Survey(name='{self.config.name}', "
            f"substances={self.n_substances}, "
            f"layers={self.config.packaging.n_layers}, "
            f"i_ref={self._i_ref})"
        )

    def save(self, path: Union[str, Path]) -> None:
        """
        Save results to NPZ file.

        Parameters
        ----------
        path : str or Path
            Output path.
        """
        if self._results is None:
            raise RuntimeError("Call compute() first")

        save_results(
            path,
            pdf_bin_centers=self._results['pdf_bin_centers'],
            pdf=self._results['pdf'],
            cdf=self._results['cdf'],
            CF_tensor=self._results['CF_tensor'],
            CF_samples=self._results['CF_samples'],
            weights=self._results['weights'],
            Fo=self._results['Fo'],
            D_vals=self._results['D_vals'],
            k_vals=self._results['k_vals'],
            k0_vals=self._results['k0_vals'],
            conc_vals=self._results['conc_vals'],
            conc_weights=self._results['conc_weights'],
            time_vals=self._results['time_vals'],
            time_weights=self._results['time_weights'],
        )

    def save_manifest(self, path: Union[str, Path]) -> None:
        """
        Save manifest (metadata) to JSON file.

        Parameters
        ----------
        path : str or Path
            Output path.
        """
        fp = fingerprint_probability(
            substances=self._substances,
            packaging=self.config.packaging,
            time_prior=self.config.time_prior,
            conc_prior=self.config.conc_prior,
            i_ref=self._i_ref,
        )

        save_manifest(
            path,
            config=self.config,
            substances=self._substances,
            fingerprint=fp,
            i_ref=self._i_ref,
            ref_layer_details=self._ref_layer_details,
            cache_stats=self._curve_cache.stats,
        )

    # =========================================================================
    # Resume Support
    # =========================================================================

    def save_state(self) -> str:
        """
        Save current state for resume capability.

        Returns
        -------
        str
            Survey ID that can be used to resume.
        """
        completed = [s.id for s in self._substances if s.id in self._curves]
        pending = [s.id for s in self._substances if s.id not in self._curves]

        state = SurveyState(
            survey_id=self._survey_id,
            config_fingerprint=fingerprint_survey_config(self.config),
            total_substances=self.n_substances,
            completed_substances=completed,
            pending_substances=pending,
        )
        self._state_manager.save(state)
        return self._survey_id

    @classmethod
    def load_state(cls, survey_id: str, config: SurveyConfig) -> Optional["Survey"]:
        """
        Load survey from saved state.

        Parameters
        ----------
        survey_id : str
            Survey ID from save_state().
        config : SurveyConfig
            Survey configuration (must match original).

        Returns
        -------
        Survey or None
            Resumed survey, or None if state not found.
        """
        state_manager = SurveyStateManager(config.cache_dir)
        state = state_manager.load(survey_id)
        if state is None:
            return None

        # Verify config matches
        current_fp = fingerprint_survey_config(config)
        if current_fp != state.config_fingerprint:
            raise ValueError("Config fingerprint mismatch. Cannot resume with different config.")

        # Create survey without substances
        survey = cls(config=config, substances=[])
        survey._survey_id = survey_id

        return survey

    # =========================================================================
    # Comparison with Fixtures
    # =========================================================================

    def compare(
        self,
        reference_path: Union[str, Path],
        tol_pdf_l1: float = 1e-3,
        tol_q_abs: float = 1e-3,
    ) -> Dict[str, Any]:
        """
        Compare results with reference fixture.

        Parameters
        ----------
        reference_path : str or Path
            Path to reference NPZ file.
        tol_pdf_l1 : float
            Tolerance for PDF L1 distance.
        tol_q_abs : float
            Tolerance for quantile absolute difference.

        Returns
        -------
        dict
            Comparison results with pass/fail status.
        """
        if self._results is None:
            raise RuntimeError("Call compute() first")

        ref = np.load(reference_path)
        ref_pdf = ref['pdf']
        ref_centers = ref['pdf_bin_centers']
        ref_cdf = ref['cdf']

        # Interpolate to common grid if needed
        our_pdf = self._results['pdf']
        our_centers = self._results['pdf_bin_centers']
        our_cdf = self._results['cdf']

        # PDF L1 distance (approximate)
        if len(our_centers) == len(ref_centers):
            dx = np.diff(our_centers).mean()
            pdf_l1 = np.sum(np.abs(our_pdf - ref_pdf)) * dx
        else:
            # Interpolate
            our_pdf_interp = np.interp(ref_centers, our_centers, our_pdf)
            dx = np.diff(ref_centers).mean()
            pdf_l1 = np.sum(np.abs(our_pdf_interp - ref_pdf)) * dx

        # Quantile comparison
        def get_quantile(cdf, centers, q):
            idx = np.searchsorted(cdf, q)
            return centers[min(idx, len(centers) - 1)]

        q50_ref = get_quantile(ref_cdf, ref_centers, 0.50)
        q95_ref = get_quantile(ref_cdf, ref_centers, 0.95)
        q99_ref = get_quantile(ref_cdf, ref_centers, 0.99)

        q50_ours = self.quantile(0.50)
        q95_ours = self.quantile(0.95)
        q99_ours = self.quantile(0.99)

        q50_diff = abs(q50_ours - q50_ref)
        q95_diff = abs(q95_ours - q95_ref)
        q99_diff = abs(q99_ours - q99_ref)

        # Pass/fail
        pdf_pass = pdf_l1 < tol_pdf_l1
        q50_pass = q50_diff < tol_q_abs
        q95_pass = q95_diff < tol_q_abs
        q99_pass = q99_diff < tol_q_abs
        all_pass = pdf_pass and q50_pass and q95_pass and q99_pass

        return {
            'pass': all_pass,
            'pdf_l1': pdf_l1,
            'pdf_l1_pass': pdf_pass,
            'q50_ref': q50_ref,
            'q50_ours': q50_ours,
            'q50_diff': q50_diff,
            'q50_pass': q50_pass,
            'q95_ref': q95_ref,
            'q95_ours': q95_ours,
            'q95_diff': q95_diff,
            'q95_pass': q95_pass,
            'q99_ref': q99_ref,
            'q99_ours': q99_ours,
            'q99_diff': q99_diff,
            'q99_pass': q99_pass,
        }
