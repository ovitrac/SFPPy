"""
survey/survey_v2.py — Bilayer and Two-Step Survey Class
=======================================================

Extends Survey for:
- Bilayer packaging (two layers, both contaminated)
- Two-step contact (chained temperatures with independent time priors)
- Combined bilayer + two-step

Falls back to the monolayer Survey for single-layer single-step scenarios.
Uses separate bilayer cache (survey/workers_v2.py + BilayerCache).

@project: SFPPy — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

import warnings
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple, Union
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict

import numpy as np

from survey.models import (
    LayerSpec, PackagingSpec, SubstanceSpec, PriorSpec,
    SurveyConfig, SubstanceModel, select_reference_layer,
)
from survey.priors import discretize_prior, compute_weight_matrix
from survey.io import load_scenario, load_substances_from_scenario
from survey.reporting import ProgressReporter
from survey.workers_v2 import (
    BilayerCurveKey, TwoStepSurfaceKey, BilayerCache,
    solve_bilayer_master_curve, solve_twostep_master_surface,
    worker_bilayer_curve, worker_twostep_surface,
)
from survey.workers_v3 import (
    BilayerCurveKeyV3, TwoStepSurfaceKeyV3, BilayerCacheV3,
    MonolayerTwoStepKeyV3,
    MultilayerCurveKeyV3, MultilayerSurfaceKeyV3,
    worker_bilayer_curve_v3, worker_twostep_surface_v3,
    worker_monolayer_twostep_surface_v3,
    worker_multilayer_curve_v3, worker_multilayer_surface_v3,
)
from survey.survey import Survey  # Inherit monolayer capabilities


class BilayerSurvey(Survey):
    """
    Survey extension for bilayer packaging and two-step contact.

    Inherits from Survey for:
    - Substance management (add/remove, inference)
    - PDF/CDF computation
    - Quantile extraction
    - Save/load

    Overrides:
    - Master curve computation (bilayer solver)
    - CF tensor computation (2D surface for two-step)
    """

    def __init__(
        self,
        config: SurveyConfig,
        substances: List[SubstanceSpec],
        substance_model: Optional[SubstanceModel] = None,
        focal_layer: int = 0,
        steps: Optional[List[Dict[str, Any]]] = None,
        step2_time_prior: Optional[PriorSpec] = None,
        use_v3: bool = False,
    ):
        """
        Parameters
        ----------
        config : SurveyConfig
            Survey configuration (must have multilayer packaging).
        substances : List[SubstanceSpec]
            Substances for the focal layer's family.
        focal_layer : int
            Which layer (0 or 1) is the source of this family's substances.
        steps : list of dict, optional
            Contact steps: [{'temperature_degC': T1}, {'temperature_degC': T2}].
        step2_time_prior : PriorSpec, optional
            Time prior for step 2 (if two-step).
        use_v3 : bool
            When True, route master-curve computation through
            survey/workers_v3.py (per-t1 chained solves,
            multi-t2 resume, equilibrium short-circuit, profile
            propagation via >>). Default False preserves the v2 baseline.
            v2 and v3 cache entries never collide (distinct key hashes
            via worker_version field).
        """
        self._focal_layer = focal_layer
        self._steps = steps
        self._step2_time_prior = step2_time_prior
        self._use_v3 = use_v3

        super().__init__(config=config, substances=substances,
                         substance_model=substance_model)

        self._bilayer_cache = (BilayerCacheV3 if use_v3 else BilayerCache)(config.cache_dir)

    @property
    def is_twostep(self) -> bool:
        return self._steps is not None and len(self._steps) >= 2

    @property
    def is_bilayer(self) -> bool:
        return self.config.packaging.n_layers >= 2

    @classmethod
    def from_scenario(cls, path: Union[str, Path],
                       use_v3: bool = False) -> "BilayerSurvey":
        """
        Create BilayerSurvey from a bilayer scenario YAML.

        `use_v3=True` routes to workers_v3 (optimised two-step workers).
        """
        import yaml
        with open(path) as f:
            raw = yaml.safe_load(f)

        config = load_scenario(path)
        substances = load_substances_from_scenario(path)

        focal_layer = raw.get('family', {}).get('focal_layer', 0)
        steps = raw.get('physics', {}).get('steps', None)

        step2_prior = None
        priors = raw.get('priors', {})
        if 'step2_time_s' in priors:
            p2 = priors['step2_time_s']
            tri = p2.get('triangular', {})
            grid = p2.get('grid', {})
            step2_prior = PriorSpec(
                mode=tri.get('mode', 180),
                max_val=tri.get('max', 360),
                n_low=grid.get('nlow', 5),
                n_high=grid.get('nhigh', 5),
            )

        obj = cls(
            config=config,
            substances=substances,
            focal_layer=focal_layer,
            steps=steps,
            step2_time_prior=step2_prior,
            use_v3=use_v3,
        )
        # Overpack-origin scenarios (converter `--overpack`, NOTE §23) self-configure the physical engine:
        # the contaminant starts on the OUTER (overpack) layer and that layer is dropped at the oven, the
        # food-contact conduit profile carried (the Δ_oven term). Primary scenarios are untouched.
        meta = raw.get('meta', {}) or {}
        if meta.get('role') == 'overpack_origin':
            obj._source_layer = int(meta.get('source_layer', -1))
            obj._overpack_removed = bool(meta.get('overpack_removed', True))
            obj._physical_engine = True
        return obj

    def _get_layer_D(self, layer_idx: int, substance: SubstanceSpec,
                     T_degC: Optional[float] = None) -> float:
        """Get diffusivity for a substance in a specific layer at temperature T."""
        from patankar.property import Dpiringer
        layer = self.config.packaging.layers[layer_idx]
        T = T_degC if T_degC is not None else layer.temperature_degC
        try:
            return Dpiringer.evaluate(layer.polymer, substance.mass_g_mol, T)
        except Exception:
            return substance.D or 1e-14

    def build_master_curves(
        self,
        parallel: bool = True,
        max_workers: Optional[int] = None,
    ) -> Dict[str, int]:
        """Build bilayer master curves (1D or 2D surface for two-step).

        A monolayer *two-step* scenario (one layer, ``physics.steps`` with
        two entries) is NOT a monolayer single-step case: it must build the
        2-D chained surface g(Fo1,Fo2) like a bilayer two-step, with a
        one-layer stack. Only a genuine monolayer 1-step falls back to the
        base Survey. Two-step is orthogonal to layer count.

        Physical-engine flag (WP-physical P2): when ``self._physical_engine`` is True, NO Fo-surface is
        built — the physical query (`_compute_cf_tensor`) solves the chain at the actual time prior.
        """
        if getattr(self, "_physical_engine", False):
            return {'hits': 0, 'misses': 0, 'physical_engine': True}

        if not self.is_bilayer and not (self.is_twostep and self._step2_time_prior):
            # Genuine monolayer single-step → base Survey.
            return super().build_master_curves(parallel=parallel,
                                                max_workers=max_workers)

        tasks = self._build_bilayer_tasks()
        if not tasks:
            return {'hits': 0, 'misses': 0}

        self._bilayer_cache._stats = {"hits": 0, "misses": 0}
        prog = ProgressReporter(total_steps=len(tasks), desc="Bilayer Curves")

        if self._use_v3:
            worker_1d = worker_bilayer_curve_v3
            worker_2d = worker_twostep_surface_v3
        else:
            worker_1d = worker_bilayer_curve
            worker_2d = worker_twostep_surface

        # Monolayer two-step uses the NATIVE one-layer max-D-clock solver
        # (not the bilayer two-step surface). Tasks were built by
        # _build_monolayer_twostep_tasks() with MonolayerTwoStepKeyV3.
        if self.is_twostep and not self.is_bilayer:
            worker_2d = worker_monolayer_twostep_surface_v3

        # n>=3 layers: tasks were built by _build_multilayer_tasks()
        # with Multilayer*KeyV3 — dispatch to the chain-engine workers.
        if self.config.packaging.n_layers >= 3:
            worker_1d = worker_multilayer_curve_v3
            worker_2d = worker_multilayer_surface_v3

        if parallel and len(tasks) > 1:
            worker_fn = worker_2d if self.is_twostep else worker_1d
            with ProcessPoolExecutor(max_workers=max_workers) as exe:
                futures = {
                    exe.submit(worker_fn, t): t["substance_id"]
                    for t in tasks
                }
                for fut in as_completed(futures):
                    result = fut.result()
                    sid = futures[fut]
                    self._store_curve_result(sid, result)
                    prog.update()
                    if result["status"] == "hit":
                        prog.hit()
                    else:
                        prog.miss()
        else:
            for task in tasks:
                if self.is_twostep:
                    result = worker_2d(task)
                else:
                    result = worker_1d(task)
                self._store_curve_result(task["substance_id"], result)
                prog.update()
                if result["status"] == "hit":
                    prog.hit()
                else:
                    prog.miss()

        prog.done()
        return self._bilayer_cache.stats

    def _store_curve_result(self, substance_id: str, result: Dict):
        """Store curve result (1D or 2D)."""
        if "surface" in result:
            self._curves[substance_id] = (
                result["fo1"], result["fo2"], result["surface"]
            )
        else:
            self._curves[substance_id] = (result["fo"], result["cf"])

    def _build_bilayer_tasks(self) -> List[Dict[str, Any]]:
        """Build task payloads for bilayer curve computation."""
        tasks = []
        pkg = self.config.packaging
        layers = pkg.layers

        if len(layers) < 2:
            if self.is_twostep and self._step2_time_prior:
                # Monolayer two-step: one-layer chained surface (not the
                # monolayer single-step fallback — that dropped step 2).
                return self._build_monolayer_twostep_tasks()
            return self._build_curve_tasks()  # genuine monolayer 1-step

        if len(layers) >= 3:
            # n>=3 layers route to the real-units chain
            # engine with tuple-valued keys. The n == 2 path below is untouched
            # (byte-identical keys → warm bilayer g-cache preserved).
            return self._build_multilayer_tasks()

        lay0 = layers[0]
        lay1 = layers[1]

        # Time prior
        time_vals, _ = discretize_prior(self.config.time_prior)
        t1_max = float(np.max(time_vals))

        # Step 2 time if two-step
        if self.is_twostep and self._step2_time_prior:
            t2_vals, _ = discretize_prior(self._step2_time_prior)
            t2_max = float(np.max(t2_vals))
            T1 = self._steps[0]['temperature_degC']
            T2 = self._steps[1]['temperature_degC']
        else:
            t2_max = None
            T1 = pkg.contact_temperature_degC
            T2 = None

        for substance in self._substances:
            D_0_T1 = self._get_layer_D(0, substance, T1)
            D_1_T1 = self._get_layer_D(1, substance, T1)

            k = substance.k or 1.0
            k0 = substance.k0 or 1.0

            # Reference layer for Fo timescale
            l_ref = lay0.thickness_m
            D_ref = D_0_T1
            perm_0 = D_0_T1 / (k * lay0.thickness_m) if k * lay0.thickness_m > 0 else float('inf')
            perm_1 = D_1_T1 / (k * lay1.thickness_m) if k * lay1.thickness_m > 0 else float('inf')
            if perm_1 < perm_0:
                l_ref = lay1.thickness_m
                D_ref = D_1_T1

            if D_ref <= 0:
                D_ref = max(D_0_T1, D_1_T1, 1e-30)

            Fo1_max = (D_ref * t1_max / l_ref**2) * self.config.fo_max_factor

            if self.is_twostep and t2_max is not None:
                D_0_T2 = self._get_layer_D(0, substance, T2)
                D_1_T2 = self._get_layer_D(1, substance, T2)
                D_ref_T2 = D_0_T2 if perm_0 <= perm_1 else D_1_T2
                if D_ref_T2 <= 0:
                    D_ref_T2 = max(D_0_T2, D_1_T2, 1e-30)
                Fo2_max = (D_ref_T2 * t2_max / l_ref**2) * self.config.fo_max_factor

                key_cls = TwoStepSurfaceKeyV3 if self._use_v3 else TwoStepSurfaceKey
                key = key_cls(
                    polymer_1=lay0.polymer, D_1_T1=D_0_T1, D_1_T2=D_0_T2,
                    k_1=k, l_1_m=lay0.thickness_m, C0_1=1.0 if self._focal_layer == 0 else 0.0,
                    polymer_2=lay1.polymer, D_2_T1=D_1_T1, D_2_T2=D_1_T2,
                    k_2=k, l_2_m=lay1.thickness_m, C0_2=1.0 if self._focal_layer == 1 else 0.0,
                    k0=k0, h=pkg.h_m_s, surface_area=pkg.surface_area_m2,
                    food_volume=pkg.food_volume_m3,
                    T1_degC=T1, T2_degC=T2, CF0=pkg.cf0,
                    Fo1_max=float(Fo1_max), Fo2_max=float(Fo2_max),
                    # Real-units chain (worker v3.2) is fast enough to afford the
                    # monolayer two-step's resolution; the old 30/15 log grid left a
                    # ~2% NON-CONSERVATIVE (under-estimate) linear-interp error on the
                    # √Fo-shaped oven curve — tightened to ~1% for the regulated path.
                    n_fo1=min(self.config.n_fo, 48),
                    n_fo2=min(self.config.n_fo, 32),
                    focal_layer=self._focal_layer,
                )
                tasks.append({
                    "cache_dir": self.config.cache_dir,
                    "key_dict": asdict(key),
                    "substance_id": substance.id,
                })
            else:
                key_cls = BilayerCurveKeyV3 if self._use_v3 else BilayerCurveKey
                key = key_cls(
                    polymer_1=lay0.polymer, D_1=D_0_T1, k_1=k,
                    l_1_m=lay0.thickness_m, C0_1=1.0 if self._focal_layer == 0 else 0.0,
                    polymer_2=lay1.polymer, D_2=D_1_T1, k_2=k,
                    l_2_m=lay1.thickness_m, C0_2=1.0 if self._focal_layer == 1 else 0.0,
                    k0=k0, h=pkg.h_m_s, surface_area=pkg.surface_area_m2,
                    food_volume=pkg.food_volume_m3,
                    contact_temperature_degC=T1, CF0=pkg.cf0,
                    Fo_max=float(Fo1_max), n_fo=self.config.n_fo,
                    focal_layer=self._focal_layer,
                )
                tasks.append({
                    "cache_dir": self.config.cache_dir,
                    "key_dict": asdict(key),
                    "substance_id": substance.id,
                })

        return tasks

    def _build_multilayer_tasks(self) -> List[Dict[str, Any]]:
        """Build n>=3-layer curve/surface tasks.

        Per-layer D is substance-specific (`_get_layer_D`); the reference layer
        is the min-permeability one (same rule as `_ref_layer_params_T1_for_
        substance`, so query-side Fo scales land on the surface's own axes).
        C0 is the focal-layer indicator, exactly as in the bilayer keys.
        Requires ``use_v3=True`` (the chain-engine solvers live in workers_v3).
        """
        if not self._use_v3:
            raise NotImplementedError(
                "n>=3 layers require use_v3=True (real-units chain solvers in "
                "survey/workers_v3.py); the v2 path is bilayer-only")

        tasks: List[Dict[str, Any]] = []
        pkg = self.config.packaging
        layers = pkg.layers
        n = len(layers)

        time_vals, _ = discretize_prior(self.config.time_prior)
        t1_max = float(np.max(time_vals))

        if self.is_twostep and self._step2_time_prior:
            t2_vals, _ = discretize_prior(self._step2_time_prior)
            t2_max = float(np.max(t2_vals))
            T1 = self._steps[0]['temperature_degC']
            T2 = self._steps[1]['temperature_degC']
        else:
            t2_max = None
            T1 = pkg.contact_temperature_degC
            T2 = None

        polymers = tuple(lay.polymer for lay in layers)
        l_ms = tuple(lay.thickness_m for lay in layers)
        C0s = tuple(1.0 if self._focal_layer == i else 0.0 for i in range(n))

        for substance in self._substances:
            k = substance.k or 1.0
            k0 = substance.k0 or 1.0
            D_T1s = tuple(self._get_layer_D(i, substance, T1) for i in range(n))

            # Reference layer: min permeability D/(k·l), first-wins on ties —
            # same rule as _ref_layer_params_T1_for_substance.
            i_ref, best_perm = 0, float('inf')
            for i in range(n):
                kl = k * l_ms[i]
                perm_i = D_T1s[i] / kl if kl > 0 else float('inf')
                if perm_i < best_perm:
                    best_perm, i_ref = perm_i, i
            l_ref = l_ms[i_ref]
            D_ref = D_T1s[i_ref]
            if D_ref <= 0:
                D_ref = max(*D_T1s, 1e-30)

            Fo1_max = (D_ref * t1_max / l_ref ** 2) * self.config.fo_max_factor

            if self.is_twostep and t2_max is not None:
                D_T2s = tuple(self._get_layer_D(i, substance, T2) for i in range(n))
                D_ref_T2 = D_T2s[i_ref]
                if D_ref_T2 <= 0:
                    D_ref_T2 = max(*D_T2s, 1e-30)
                Fo2_max = (D_ref_T2 * t2_max / l_ref ** 2) * self.config.fo_max_factor
                key = MultilayerSurfaceKeyV3(
                    polymers=polymers, D_T1s=D_T1s, D_T2s=D_T2s,
                    ks=tuple(k for _ in range(n)), l_ms=l_ms, C0s=C0s,
                    k0=k0, h=pkg.h_m_s, surface_area=pkg.surface_area_m2,
                    food_volume=pkg.food_volume_m3,
                    T1_degC=T1, T2_degC=T2, CF0=pkg.cf0,
                    Fo1_max=float(Fo1_max), Fo2_max=float(Fo2_max),
                    n_fo1=min(self.config.n_fo, 48), n_fo2=min(self.config.n_fo, 32),
                    focal_layer=self._focal_layer,
                )
            else:
                key = MultilayerCurveKeyV3(
                    polymers=polymers, D_T1s=D_T1s,
                    ks=tuple(k for _ in range(n)), l_ms=l_ms, C0s=C0s,
                    k0=k0, h=pkg.h_m_s, surface_area=pkg.surface_area_m2,
                    food_volume=pkg.food_volume_m3,
                    contact_temperature_degC=T1, CF0=pkg.cf0,
                    Fo_max=float(Fo1_max), n_fo=self.config.n_fo,
                    focal_layer=self._focal_layer,
                )
            tasks.append({
                "cache_dir": self.config.cache_dir,
                "key_dict": asdict(key),
                "substance_id": substance.id,
            })

        return tasks

    def _build_monolayer_twostep_tasks(self) -> List[Dict[str, Any]]:
        """Build NATIVE monolayer two-step surface tasks (max-D common clock).

        Honours the scenario as set: one focal layer, two temperatures, two
        chained time priors → one cached surface g(Fo1,Fo2) per physical key.
        This is a *real one-layer* task (``MonolayerTwoStepKeyV3``), NOT a
        degenerate bilayer. Both Fo axes use the single common clock
        ``D_ref = max(D@T1, D@T2)``, so the normalised diffusivity is ``<= 1``
        in both steps (no stiffness; the slow storage step is a large formal
        clock value with tiny ``D_norm``, not an inflated fast horizon). The
        distinct key hash space keeps these entries from ever colliding with
        shipped bilayer surfaces.

        Requires ``use_v3=True`` (the native monolayer solver lives in
        workers_v3); the v2 path is bilayer-only.
        """
        if not self._use_v3:
            raise NotImplementedError(
                "monolayer two-step requires use_v3=True (native max-D-clock "
                "solver lives in survey/workers_v3.py); the v2 two-step path "
                "is bilayer-only")

        tasks: List[Dict[str, Any]] = []
        pkg = self.config.packaging
        lay0 = pkg.layers[0]

        time_vals, _ = discretize_prior(self.config.time_prior)
        t1_max = float(np.max(time_vals))
        t2_vals, _ = discretize_prior(self._step2_time_prior)
        t2_max = float(np.max(t2_vals))
        T1 = self._steps[0]['temperature_degC']
        T2 = self._steps[1]['temperature_degC']

        l_ref = lay0.thickness_m

        for substance in self._substances:
            D_T1 = self._get_layer_D(0, substance, T1)
            D_T2 = self._get_layer_D(0, substance, T2)
            k = substance.k or 1.0
            k0 = substance.k0 or 1.0

            # Common clock anchored to the HIGHEST D in the sequence → both
            # Fo axes on D_ref; the slow step gets D_norm<<1 (non-stiff), the
            # fast step D_norm==1 (no ×ratio horizon inflation).
            D_ref = max(D_T1, D_T2)
            if D_ref <= 0:
                D_ref = max(D_T1, D_T2, 1e-30)
            Fo1_max = (D_ref * t1_max / l_ref ** 2) * self.config.fo_max_factor
            Fo2_max = (D_ref * t2_max / l_ref ** 2) * self.config.fo_max_factor

            key = MonolayerTwoStepKeyV3(
                polymer=lay0.polymer, D_T1=D_T1, D_T2=D_T2,
                k=k, l_m=lay0.thickness_m, C0=1.0,
                k0=k0, h=pkg.h_m_s, surface_area=pkg.surface_area_m2,
                food_volume=pkg.food_volume_m3,
                T1_degC=T1, T2_degC=T2, CF0=pkg.cf0,
                Fo1_max=float(Fo1_max), Fo2_max=float(Fo2_max),
                # √Fo-spaced monolayer grids are cheap (one fast layer) → afford
                # higher resolution than the bilayer surface for ~1% accuracy.
                n_fo1=min(self.config.n_fo, 64),
                n_fo2=min(self.config.n_fo, 32),
            )
            tasks.append({
                "cache_dir": self.config.cache_dir,
                "key_dict": asdict(key),
                "substance_id": substance.id,
            })
        return tasks

    # ------------------------------------------------------------------
    # Patch P3 (2026-05-12): per-substance reference layer parameters
    # for the fo→t scale used to query the bilayer cache.
    #
    # The bilayer master surface is indexed in dimensionless Fourier
    # coordinates anchored to the min-permeability layer's diffusivity
    # at the contact temperature (`D_ref_T1`) and that same layer's
    # thickness (`l_ref`).  The full scale `α = D_ref / l_ref²` MUST
    # come from the same `i_ref` returned by the per-substance
    # permeability comparison — never D from one layer with l from
    # another.
    #
    # Before P3, `_compute_cf_tensor_*` used `substance.D / l_ref²`
    # with `l_ref` taken from `self.ref_layer` (scenario-wide), which
    # disagreed with the cache's per-substance reference layer
    # whenever `substance.D ≠ D_ref` or whenever the per-substance
    # `i_ref` differed from `self._i_ref` (fixed in the v3 survey-engine
    # consolidation).
    # ------------------------------------------------------------------

    def _ref_layer_params_T1_for_substance(
        self, substance: "SubstanceSpec"
    ) -> Tuple[int, float, float]:
        """
        Return the paired triple ``(i_ref, D_ref_T1, l_ref)`` for
        ``substance`` using the same min-permeability rule as
        ``_build_bilayer_tasks`` (`survey_v2.py:262–271`).

        ``D_ref_T1`` and ``l_ref`` are guaranteed to come from the same
        ``i_ref`` — never select one without the other.
        """
        pkg = self.config.packaging
        if len(pkg.layers) < 2:
            # Monolayer fallback — caller should not reach here in
            # bilayer paths, but be defensive.
            l_ref = pkg.layers[0].thickness_m
            T1 = pkg.contact_temperature_degC
            D_ref = self._get_layer_D(0, substance, T1)
            return 0, D_ref, l_ref

        # n-layer min-permeability loop. Strict '<' keeps
        # the earlier (food-contact-most) layer on ties — for n == 2 this is
        # BIT-IDENTICAL to the previous pairwise comparison, so all cached
        # bilayer Fo scales are unchanged.
        T1 = pkg.contact_temperature_degC
        k = substance.k or 1.0
        i_ref, D_ref, l_ref = 0, 0.0, pkg.layers[0].thickness_m
        best_perm = float('inf')
        for i, lay in enumerate(pkg.layers):
            D_i = self._get_layer_D(i, substance, T1)
            kl = k * lay.thickness_m
            perm_i = D_i / kl if kl > 0 else float('inf')
            if perm_i < best_perm:
                best_perm = perm_i
                i_ref, D_ref, l_ref = i, D_i, lay.thickness_m
        if best_perm == float('inf'):
            # All layers degenerate — keep layer 0 with its D (defensive).
            D_ref = self._get_layer_D(0, substance, T1)
        return i_ref, D_ref, l_ref

    def _is_monolayer_twostep(self) -> bool:
        return (self.is_twostep and not self.is_bilayer
                and self._steps is not None and len(self._steps) >= 2)

    def _mono_twostep_clock_scale(self, substance: "SubstanceSpec") -> float:
        """Common-clock Fo scale ``D_ref / l²`` for a monolayer two-step, with
        ``D_ref = max(D@T1, D@T2)`` — identical to ``MonolayerTwoStepKeyV3`` so
        the interpolation lands on the surface's own axes. BOTH Fo axes share
        this single scale (the max-D common clock)."""
        pkg = self.config.packaging
        l_ref = pkg.layers[0].thickness_m
        T1 = self._steps[0].get('temperature_degC')
        T2 = self._steps[1].get('temperature_degC')
        D_T1 = self._get_layer_D(0, substance, T1)
        D_T2 = self._get_layer_D(0, substance, T2)
        D_ref = max(D_T1, D_T2)
        if D_ref <= 0:
            D_ref = max(D_T1, D_T2, 1e-30)
        return D_ref / (l_ref ** 2)

    def _ref_fo_scale_T1_for_substance(
        self, substance: "SubstanceSpec"
    ) -> float:
        """High-level: returns ``α_T1 = D_ref_T1 / l_ref²`` from the
        same ``i_ref`` as the cache build. Monolayer two-step uses the
        max-D common clock (both axes share ``max(D@T1,D@T2)/l²``)."""
        if self._is_monolayer_twostep():
            return self._mono_twostep_clock_scale(substance)
        _, D_ref_T1, l_ref = self._ref_layer_params_T1_for_substance(substance)
        return D_ref_T1 / (l_ref ** 2)

    def _ref_fo_scale_T2_for_substance(
        self, substance: "SubstanceSpec"
    ) -> float:
        """High-level: returns ``α_T2 = D_ref_T2 / l_ref²`` using the
        SAME ``i_ref`` as ``_ref_layer_params_T1_for_substance`` and
        the step-2 temperature from ``self._steps[1]``. Monolayer two-step
        uses the max-D common clock (same scale as T1)."""
        if self._is_monolayer_twostep():
            return self._mono_twostep_clock_scale(substance)
        i_ref, _, l_ref = self._ref_layer_params_T1_for_substance(substance)
        T2 = None
        if self._steps and len(self._steps) >= 2:
            T2 = self._steps[1].get('temperature_degC')
        D_ref_T2 = self._get_layer_D(i_ref, substance, T2)
        return D_ref_T2 / (l_ref ** 2)

    def _compute_cf_tensor(self) -> Dict[str, np.ndarray]:
        """
        Compute CF tensor for bilayer scenarios.

        For single-step: CF(t, CP0) = g(Fo(t)) × CP0  (same as monolayer)
        For two-step: CF(t1, t2, CP0) = g(Fo1(t1), Fo2(t2)) × CP0  (3D tensor)

        Physical-engine flag (WP-physical P2): when ``self._physical_engine`` is True, the family CF
        tensor is computed by `survey.physical_query.physical_cf_tensor` — the real-units chain
        evaluated at the ACTUAL discretised time prior (no Fo-surface, no interpolation). Default OFF;
        the legacy Fo-surface path below is unchanged. Validated full-key-parity drop-in
        (T11: vs legacy mono-1 1.7e-7, bi-1 2.9e-16, bi-2 0.43%). Covers all cases (mono-1 = N=1).
        """
        if getattr(self, "_physical_engine", False):
            from survey.physical_query import physical_cf_tensor
            return physical_cf_tensor(self)

        if not self.is_bilayer and not (self.is_twostep and self._step2_time_prior):
            # Genuine monolayer single-step → base Survey curve path.
            return super()._compute_cf_tensor()

        # Discretize priors
        time_vals, time_w = discretize_prior(self.config.time_prior)
        conc_vals, conc_w = discretize_prior(self.config.conc_prior)

        # Reference layer for Fo mapping (monolayer two-step: the only layer)
        layer = self.ref_layer
        l2 = layer.thickness_m ** 2

        if self.is_twostep and self._step2_time_prior:
            return self._compute_cf_tensor_twostep(
                time_vals, time_w, conc_vals, conc_w, l2)
        else:
            return self._compute_cf_tensor_bilayer_1step(
                time_vals, time_w, conc_vals, conc_w, l2)

    def per_substance_cf_tensors(self, step: str = "full") -> Dict[str, Dict[str, np.ndarray]]:
        """Bilayer per-substance decomposition; ``Σ_i CF_i == family``.

        Covers 1-step (2-D ``outer(g_i, conc)``) and two-step (3-D
        ``einsum('ij,k', g_i_surface, conc)``), reusing the same per-substance
        Fo scales and master curves/surfaces as ``_compute_cf_tensor*``. The
        contribution factor (p_i / w_norm_i) is applied after the curve, so the
        concentration enters post-simulation exactly as in the family path.

        ``step`` selects the exposure chain (storage/full step split):
          - ``"full"``    → step 1+2, storage then oven (default; bit-identical
            to the historical behaviour).
          - ``"storage"`` → step 1 alone, the surface sliced at Fo2=0 (the same
            g1 construction as ``step_resolved_cf_tensors``). For a genuine
            1-step scenario storage ≡ full (no oven), so ``step`` is a no-op there.

        Physical-engine flag (WP-physical P2): routed through the real-units chain when set.
        """
        if step not in ("full", "storage"):
            raise ValueError(f"step must be 'full' or 'storage', got {step!r}")
        if getattr(self, "_physical_engine", False):
            from survey.physical_query import physical_per_substance_cf_tensors
            return physical_per_substance_cf_tensors(self, step=step)
        # A two-step MONOLAYER is not a bilayer but its master curve is a 2-D (fo1,fo2,surface), so it
        # must use the two-step surface branch below — NOT the base 1-D path (which unpacks a 2-tuple and
        # raises on the 3-tuple surface). Only a genuine monolayer 1-step falls back to the base Survey.
        if not self.is_bilayer and not (self.is_twostep and self._step2_time_prior):
            return super().per_substance_cf_tensors()  # 1-step: storage ≡ full
        time_vals, time_w = discretize_prior(self.config.time_prior)
        conc_vals, conc_w = discretize_prior(self.config.conc_prior)
        factor = self._per_substance_factors()
        subs = self._substances
        out: Dict[str, Dict[str, np.ndarray]] = {}

        def _store(cas, cf, extra):
            if cas in out:
                out[cas]['CF_tensor'] = out[cas]['CF_tensor'] + cf
                out[cas]['CF_samples'] = out[cas]['CF_tensor'].ravel()
            else:
                d = {'CF_tensor': cf, 'CF_samples': cf.ravel()}
                d.update(extra)
                out[cas] = d

        if self.is_twostep and self._step2_time_prior:
            time2_vals, time2_w = discretize_prior(self._step2_time_prior)
            scale_T1 = [self._ref_fo_scale_T1_for_substance(s) for s in subs]
            scale_T2 = [self._ref_fo_scale_T2_for_substance(s) for s in subs]
            from scipy.interpolate import RegularGridInterpolator
            n_t1, n_t2 = len(time_vals), len(time2_vals)
            if step == "storage":
                # Step 1 alone: slice each substance surface at Fo2=0 (identical g1
                # construction to step_resolved_cf_tensors), giving a single-step-format
                # 2-D (n_t1, n_cp0) tensor so the downstream combiners treat it as 1-step.
                wmat = compute_weight_matrix(time_w, conc_w).ravel()
                extra_s = {'weights': wmat, 'time_vals': time_vals, 'time_weights': time_w,
                           'conc_vals': conc_vals, 'conc_weights': conc_w}
                for j, s in enumerate(subs):
                    fo1_grid, fo2_grid, surface = self._curves[s.id][:3]
                    interp = RegularGridInterpolator(
                        (fo1_grid, fo2_grid), surface, method='linear',
                        bounds_error=False, fill_value=None)
                    Fo1 = np.clip(time_vals * scale_T1[j], fo1_grid[0], fo1_grid[-1])
                    fo2_0 = max(0.0, float(fo2_grid[0]))
                    g1 = interp(np.column_stack([Fo1, np.full(n_t1, fo2_0)]))
                    cf = np.outer(factor[j] * g1, conc_vals)
                    cas = getattr(s, 'cas', None) or getattr(s, 'id', f"idx{j}")
                    _store(cas, cf, dict(extra_s,
                                         exchangeable=bool(getattr(s, 'exchangeable', True)),
                                         weight=float(getattr(s, 'weight', 1.0))))
                return out
            w3 = np.einsum('i,j,k->ijk', time_w, time2_w, conc_w)
            w3 /= w3.sum()
            # step2_active=True: the 3-D (n_t1, n_t2, n_cp0) per-substance tensor must
            # be marginalised over the oven-time (t2) prior by combine_sources_shared_t;
            # without it the aggregator iterates t1 only and raises "t2_idx required
            # for step-2-active source with 3-D CF_tensor" (mirrors the physical path).
            extra = {'weights': w3.ravel(), 'time_vals': time_vals,
                     'time_weights': time_w, 'time2_vals': time2_vals,
                     'time2_weights': time2_w, 'conc_vals': conc_vals,
                     'conc_weights': conc_w, 'step2_active': True}
            for j, s in enumerate(subs):
                fo1_grid, fo2_grid, surface = self._curves[s.id][:3]
                Fo1 = np.clip(time_vals * scale_T1[j], fo1_grid[0], fo1_grid[-1])
                Fo2 = np.clip(time2_vals * scale_T2[j], fo2_grid[0], fo2_grid[-1])
                interp = RegularGridInterpolator(
                    (fo1_grid, fo2_grid), surface, method='linear',
                    bounds_error=False, fill_value=None)
                M1, M2 = np.meshgrid(Fo1, Fo2, indexing='ij')
                g2d = interp(np.column_stack([M1.ravel(), M2.ravel()])).reshape(n_t1, n_t2)
                cf = np.einsum('ij,k->ijk', factor[j] * g2d, conc_vals)
                cas = getattr(s, 'cas', None) or getattr(s, 'id', f"idx{j}")
                _store(cas, cf, dict(extra,
                                     exchangeable=bool(getattr(s, 'exchangeable', True)),
                                     weight=float(getattr(s, 'weight', 1.0))))
        else:
            scale_T1 = [self._ref_fo_scale_T1_for_substance(s) for s in subs]
            Fo = np.outer(time_vals, np.array(scale_T1, dtype=float))
            wmat = compute_weight_matrix(time_w, conc_w).ravel()
            extra = {'weights': wmat, 'time_vals': time_vals, 'time_weights': time_w,
                     'conc_vals': conc_vals, 'conc_weights': conc_w}
            for j, s in enumerate(subs):
                fo_grid, cf_over_cp0 = self._curves[s.id][0], self._curves[s.id][1]
                g = np.interp(np.clip(Fo[:, j], 0.0, None), fo_grid, cf_over_cp0)
                cf = np.outer(factor[j] * g, conc_vals)
                cas = getattr(s, 'cas', None) or getattr(s, 'id', f"idx{j}")
                _store(cas, cf, dict(extra,
                                     exchangeable=bool(getattr(s, 'exchangeable', True)),
                                     weight=float(getattr(s, 'weight', 1.0))))
        return out

    def step_resolved_cf_tensors(self) -> Dict[str, np.ndarray]:
        """Step-resolved family CF tensors for a two-step scenario, all sliced
        from the cached per-substance master surfaces (no re-sim):

          CF1  = storage alone       = surface at Fo2=0  → (n_t1, n_c)
          CF2  = oven heating alone  = surface at Fo1=0  → (n_t2, n_c)
          CF12 = storage then oven   = full surface      → (n_t1, n_t2, n_c)

        Heating effect: approx = CF12 - CF1 (broadcast over t2),
                        exact  = CF12 - CF2 (broadcast over t1); both paired
        per cell, then quantiled with w12.
        """
        if not (self.is_twostep and self._step2_time_prior):
            raise ValueError("step_resolved_cf_tensors requires a two-step scenario")
        if getattr(self, "_physical_engine", False):
            from survey.physical_query import physical_step_resolved_cf_tensors
            return physical_step_resolved_cf_tensors(self)
        from scipy.interpolate import RegularGridInterpolator
        t1_vals, t1_w = discretize_prior(self.config.time_prior)
        t2_vals, t2_w = discretize_prior(self._step2_time_prior)
        conc_vals, conc_w = discretize_prior(self.config.conc_prior)
        subs = self._substances
        scale_T1 = [self._ref_fo_scale_T1_for_substance(s) for s in subs]
        scale_T2 = [self._ref_fo_scale_T2_for_substance(s) for s in subs]
        n_t1, n_t2 = len(t1_vals), len(t2_vals)
        surf_full, g1, g2 = [], [], []
        for j, s in enumerate(subs):
            fo1_grid, fo2_grid, surface = self._curves[s.id][:3]
            interp = RegularGridInterpolator((fo1_grid, fo2_grid), surface,
                                             method='linear', bounds_error=False,
                                             fill_value=None)
            Fo1 = np.clip(t1_vals * scale_T1[j], fo1_grid[0], fo1_grid[-1])
            Fo2 = np.clip(t2_vals * scale_T2[j], fo2_grid[0], fo2_grid[-1])
            fo2_0 = max(0.0, float(fo2_grid[0]))   # storage alone → no oven
            fo1_0 = max(0.0, float(fo1_grid[0]))   # oven alone → no storage
            M1, M2 = np.meshgrid(Fo1, Fo2, indexing='ij')
            surf_full.append(interp(np.column_stack([M1.ravel(), M2.ravel()])).reshape(n_t1, n_t2))
            g1.append(interp(np.column_stack([Fo1, np.full(n_t1, fo2_0)])))
            g2.append(interp(np.column_stack([np.full(n_t2, fo1_0), Fo2])))
        g12 = self._combine_substance_surfaces(surf_full)            # (n_t1, n_t2)
        g1c = self._combine_substance_curves(np.asarray(g1).T)        # (n_t1,)
        g2c = self._combine_substance_curves(np.asarray(g2).T)        # (n_t2,)
        w12 = np.einsum('i,j,k->ijk', t1_w, t2_w, conc_w); w12 /= w12.sum()
        return {
            'CF1': np.outer(g1c, conc_vals), 'CF2': np.outer(g2c, conc_vals),
            'CF12': np.einsum('ij,k->ijk', g12, conc_vals),
            'w1': compute_weight_matrix(t1_w, conc_w).ravel(),
            'w2': compute_weight_matrix(t2_w, conc_w).ravel(),
            'w12': w12.ravel(),
            't1_vals': t1_vals, 't2_vals': t2_vals, 'conc_vals': conc_vals,
        }

    def _compute_cf_tensor_bilayer_1step(
        self, time_vals, time_w, conc_vals, conc_w, l2
    ) -> Dict[str, np.ndarray]:
        """CF tensor for bilayer single-step (same structure as monolayer)."""
        n_sub = len(self._substances)
        n_t = len(time_vals)

        # Patch P3 (2026-05-12): use the per-substance fo→t scale
        # `α = D_ref / l_ref²` from the same `i_ref` as the cache build.
        # The previous mapping `substance.D / l_ref²` (with scenario-
        # wide `l_ref`) disagreed with the cached surface coordinates.
        # (fixed in the v3 survey-engine consolidation)
        scale_T1 = np.array(
            [self._ref_fo_scale_T1_for_substance(s) for s in self._substances],
            dtype=float,
        )
        # `D_vals` is reported alongside the tensor for back-compat
        # consumers; the actual scale used for interp is `scale_T1`.
        D_vals = np.array([s.D or 1e-14 for s in self._substances], dtype=float)
        Fo = np.outer(time_vals, scale_T1)

        # Interpolate master curves
        g_tj = np.empty((n_t, n_sub), dtype=float)
        for j, substance in enumerate(self._substances):
            curve = self._curves[substance.id]
            fo_grid, cf_over_cp0 = curve[0], curve[1]
            g_tj[:, j] = np.interp(
                np.clip(Fo[:, j], 0.0, None), fo_grid, cf_over_cp0)

        # Exchangeable/non-exchangeable weighting (same logic as Survey)
        g_combined = self._combine_substance_curves(g_tj)

        # CF = g_combined(t) × CP0
        cf_family = np.outer(g_combined, conc_vals)

        # Family-level mass-conservation check (proper probabilistic weights).
        # Per-substance conservation is enforced in the solver; here we bound the
        # COMBINED family curve by the same weighting applied to the per-substance
        # ceiling V_focal/V_F (focal layer = the only loaded layer). Uses V_focal,
        # NOT the total bilayer volume, per the kernel mass rule.
        self._check_family_mass_conservation(g_combined)

        substance_weights = np.array([
            getattr(s, 'weight', 1.0) for s in self._substances], dtype=float)
        w_sum = substance_weights.sum()
        if w_sum > 0:
            substance_weights /= w_sum

        return {
            'CF_tensor': cf_family,
            'CF_samples': cf_family.ravel(),
            'weights': compute_weight_matrix(time_w, conc_w).ravel(),
            'Fo': Fo,
            'D_vals': D_vals,
            'conc_vals': conc_vals,
            'conc_weights': conc_w,
            'time_vals': time_vals,
            'time_weights': time_w,
            'substance_weights': substance_weights,
            'n_exchangeable': sum(1 for s in self._substances
                                  if getattr(s, 'exchangeable', True)),
            'n_non_exchangeable': sum(1 for s in self._substances
                                      if not getattr(s, 'exchangeable', True)),
        }

    def _compute_cf_tensor_twostep(
        self, time1_vals, time1_w, conc_vals, conc_w, l2
    ) -> Dict[str, np.ndarray]:
        """
        CF tensor for two-step: 3D (t1 × t2 × CP0).

        g(Fo1, Fo2) is a 2D surface per substance.
        The combined g is the weighted average (exchangeable) or
        weighted sum (non-exchangeable) of per-substance surfaces.
        """
        time2_vals, time2_w = discretize_prior(self._step2_time_prior)

        # Patch P3 (2026-05-12): per-substance fo→t scale via
        # `(i_ref, D_ref, l_ref)` triple — see notes on
        # `_ref_layer_params_T1_for_substance`.
        scale_T1 = np.array(
            [self._ref_fo_scale_T1_for_substance(s) for s in self._substances],
            dtype=float,
        )
        scale_T2 = np.array(
            [self._ref_fo_scale_T2_for_substance(s) for s in self._substances],
            dtype=float,
        )
        # `D_vals` reported alongside the tensor for back-compat
        # consumers; the actual scales used for interp are `scale_T*`.
        D_vals = np.array([s.D or 1e-14 for s in self._substances], dtype=float)
        n_sub = len(self._substances)
        n_t1 = len(time1_vals)
        n_t2 = len(time2_vals)

        # For each substance, interpolate the 2D surface
        g_surfaces = []  # List of (n_t1, n_t2) arrays
        for j, substance in enumerate(self._substances):
            curve = self._curves[substance.id]
            fo1_grid, fo2_grid, surface = curve[0], curve[1], curve[2]

            # Map real times to Fo using the per-substance scale (P3).
            Fo1 = time1_vals * scale_T1[j]
            Fo2 = time2_vals * scale_T2[j]

            # 2D interpolation on the surface
            from scipy.interpolate import RegularGridInterpolator
            interp = RegularGridInterpolator(
                (fo1_grid, fo2_grid), surface,
                method='linear', bounds_error=False, fill_value=None,
            )
            # Build meshgrid of (Fo1, Fo2) points
            Fo1_mesh, Fo2_mesh = np.meshgrid(
                np.clip(Fo1, fo1_grid[0], fo1_grid[-1]),
                np.clip(Fo2, fo2_grid[0], fo2_grid[-1]),
                indexing='ij',
            )
            points = np.column_stack([Fo1_mesh.ravel(), Fo2_mesh.ravel()])
            g_2d = interp(points).reshape(n_t1, n_t2)
            g_surfaces.append(g_2d)

        # Combine substances (exchangeable weighted average)
        g_combined_2d = self._combine_substance_surfaces(g_surfaces)

        # Family-level mass-conservation check (same proper-weight bound as the
        # 1-step path; surface max is the relevant peak for two-step).
        self._check_family_mass_conservation(g_combined_2d)

        # CF tensor: 3D = g(t1, t2) × CP0
        # Shape: (n_t1, n_t2, n_c)
        cf_3d = np.einsum('ij,k->ijk', g_combined_2d, conc_vals)

        # Joint weights: 3D = w_t1 × w_t2 × w_c
        w_3d = np.einsum('i,j,k->ijk', time1_w, time2_w, conc_w)
        w_3d /= w_3d.sum()

        return {
            'CF_tensor': cf_3d,
            'CF_samples': cf_3d.ravel(),
            'weights': w_3d.ravel(),
            'D_vals': D_vals,
            'conc_vals': conc_vals,
            'conc_weights': conc_w,
            'time_vals': time1_vals,
            'time_weights': time1_w,
            'time2_vals': time2_vals,
            'time2_weights': time2_w,
            'substance_weights': np.ones(n_sub) / max(n_sub, 1),
            'n_exchangeable': sum(1 for s in self._substances
                                  if getattr(s, 'exchangeable', True)),
            'n_non_exchangeable': sum(1 for s in self._substances
                                      if not getattr(s, 'exchangeable', True)),
        }

    def _check_family_mass_conservation(self, g_combined) -> None:
        """Family-level mass-conservation check for the bilayer paths.

        Bounds the combined family ``CF/CP0`` (1-D curve or 2-D surface) by the
        per-substance ceiling ``V_focal/V_F`` (focal = the only loaded layer)
        propagated through the exact family-combination weighting
        (``_family_mass_ceiling``). Warns loudly on breach; the per-substance
        kernel guard in the solver is the independent first check.
        """
        from survey.workers import MASS_CONSERVATION_TOL
        pkg = self.config.packaging
        V_F = pkg.food_volume_m3
        if not (V_F > 0):
            return
        focal = getattr(self, '_focal_layer', 0) or 0
        layers = pkg.layers
        l_focal = layers[focal].thickness_m if 0 <= focal < len(layers) \
            else layers[0].thickness_m
        per_sub_ceiling = l_focal * pkg.surface_area_m2 / V_F
        family_ceiling = self._family_mass_ceiling(per_sub_ceiling)
        g_max = float(np.nanmax(g_combined))
        if g_max > family_ceiling * (1.0 + MASS_CONSERVATION_TOL):
            warnings.warn(
                f"family CF/CP0 max={g_max:.4g} exceeds weighted family mass "
                f"bound={family_ceiling:.4g} (per-substance V_focal/V_F="
                f"{per_sub_ceiling:.4g}); check substance combination weights",
                RuntimeWarning,
            )

    def _combine_substance_curves(self, g_tj: np.ndarray) -> np.ndarray:
        """
        Combine per-substance master curves into family curve.

        Same exchangeable/non-exchangeable logic as Survey._compute_cf_tensor(),
        but extracted into a reusable method.
        """
        n_t = g_tj.shape[0]
        exch_idx = [i for i, s in enumerate(self._substances)
                    if getattr(s, 'exchangeable', True)]
        non_exch_idx = [i for i, s in enumerate(self._substances)
                        if not getattr(s, 'exchangeable', True)]

        g_combined = np.zeros(n_t, dtype=float)
        n_exch = len(exch_idx)

        if n_exch > 0:
            w_exch = np.array([getattr(self._substances[i], 'weight', 1.0)
                               for i in exch_idx], dtype=float)
            w_sum = w_exch.sum()
            p_exch = w_exch / w_sum if w_sum > 0 else np.ones(n_exch) / n_exch
            g_exch = np.dot(g_tj[:, exch_idx], p_exch)
            g_combined += g_exch  # weighted average only; no 1/E (F8, 2026-06-12)

        if non_exch_idx:
            w_non = np.array([getattr(self._substances[i], 'weight', 1.0)
                              for i in non_exch_idx], dtype=float)
            w_max = w_non.max()
            w_norm = w_non / w_max if w_max > 0 else np.ones(len(non_exch_idx))
            g_non = np.dot(g_tj[:, non_exch_idx], w_norm)
            g_combined += g_non

        return g_combined

    def _combine_substance_surfaces(self, surfaces: List[np.ndarray]) -> np.ndarray:
        """
        Combine per-substance 2D surfaces into family surface.

        Same logic as _combine_substance_curves but for 2D arrays.
        """
        shape = surfaces[0].shape
        exch_idx = [i for i, s in enumerate(self._substances)
                    if getattr(s, 'exchangeable', True)]
        non_exch_idx = [i for i, s in enumerate(self._substances)
                        if not getattr(s, 'exchangeable', True)]

        g_combined = np.zeros(shape, dtype=float)
        n_exch = len(exch_idx)

        if n_exch > 0:
            w_exch = np.array([getattr(self._substances[i], 'weight', 1.0)
                               for i in exch_idx], dtype=float)
            w_sum = w_exch.sum()
            p_exch = w_exch / w_sum if w_sum > 0 else np.ones(n_exch) / n_exch
            g_exch = sum(surfaces[i] * p_exch[j]
                         for j, i in enumerate(exch_idx))
            g_combined += g_exch  # weighted average only; no 1/E (F8, 2026-06-12)

        if non_exch_idx:
            w_non = np.array([getattr(self._substances[i], 'weight', 1.0)
                              for i in non_exch_idx], dtype=float)
            w_max = w_non.max()
            w_norm = w_non / w_max if w_max > 0 else np.ones(len(non_exch_idx))
            g_non = sum(surfaces[i] * w_norm[j]
                        for j, i in enumerate(non_exch_idx))
            g_combined += g_non

        return g_combined
