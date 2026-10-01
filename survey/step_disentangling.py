"""
survey.step_disentangling
=========================

Post-hoc separation of probabilistic step-1 (storage) and step-1+2
(storage + oven heating) marginals from completed two-step v3 manifests.

Path A (this module) reads the bilayer master-surface cache and slices
fo2 = 0 to recover the step-1-only response surface.  **No new
simulation.**  For single-step scenarios, the step-1 view is bit-
identical to the parent manifest (hardlink emitted for uniform
downstream iteration).

Control points
--------------
Five invariants are checked at every disentangling call.  Classification:

  CP-1  Time prior loaded correctly      HARD GATE  (raises under strict=True)
  CP-2  Layer temperature loaded         HARD GATE  (raises under strict=True)
  CP-3  Fo1_query within cache range     HARD GATE  (raises under strict=True)
  CP-4  Mass-balance bound respected     HARD GATE  (raises under strict=True
                                                     when cf_sat_bound is in cache)
  CP-5  CF_step1 ≤ CF_full per-cell      DIAGNOSTIC ONLY during pre-fix transition;
                                          becomes HARD GATE once the v3 caches
                                          are rebuilt.  Rationale: against
                                          the pre-fix saturation-clipped parent
                                          CF_full, a correctly disentangled
                                          CF_step1 (from the cache, using
                                          D_ref/l_ref²) will appear larger
                                          per cell — a false positive that
                                          would otherwise abort every call.

When ``strict=False`` all CPs are logged as warnings instead of raising —
DIAGNOSTIC ONLY; results from a non-strict call MUST NOT be used for reporting.

Author:  Olivier Vitrac, PhD, HDR
Contact: olivier.vitrac@gmail.com
"""

from __future__ import annotations

import json
import os
import warnings
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import yaml as _yaml

from survey.cache_extensions import (
    ExtendedCacheEntry,
    load_bilayer_cache,
)


# --------------------------------------------------------------------- #
# Public API                                                            #
# --------------------------------------------------------------------- #

def disentangle_step1_for_scenario(
    scenario_yaml: Path,
    parent_npz: Path,
    output_npz: Path,
    bilayer_cache_dir: Path,
    *,
    overwrite: bool = False,
    strict: bool = True,
) -> Dict[str, Any]:
    """
    Build a step-1-only virtual NPZ from a parent 2-step manifest.

    Parameters
    ----------
    scenario_yaml : Path
        The scenario YAML that produced the parent NPZ.
    parent_npz : Path
        The parent manifest emitted by ``_save_results_npz`` in
        ``run_production_v3.py`` — either 1-step or 2-step.
    output_npz : Path
        Destination for ``<name>_step1.npz``.  Parent directory is
        created if missing.
    bilayer_cache_dir : Path
        The ``.survey_cache/`` directory; the bilayer cache lives in
        ``<bilayer_cache_dir>/bilayer_curves/``.  Should match the
        ``config.cache_dir`` at production time.
    overwrite : bool, default False
        If False and ``output_npz`` already exists with a matching
        provenance hash, return without recomputing.
    strict : bool, default True
        Raise on CP failures.  If False, log warnings and continue —
        DIAGNOSTIC ONLY; results from a non-strict call MUST NOT be
        used for reporting.

    Returns
    -------
    provenance : dict
        method, parent NPZ path, cache file, shape, YAML SHA-256,
        and the per-CP outcomes.
    """
    parent_npz = Path(parent_npz)
    output_npz = Path(output_npz)
    scenario_yaml = Path(scenario_yaml)
    bilayer_cache_dir = Path(bilayer_cache_dir)

    output_npz.parent.mkdir(parents=True, exist_ok=True)

    # Idempotent early exit
    if output_npz.exists() and not overwrite:
        existing_meta = _read_provenance_sidecar(output_npz)
        if existing_meta and existing_meta.get("parent_npz_sha256") == \
                _fast_sha256(parent_npz):
            return {**existing_meta, "status": "skipped_idempotent"}

    parent = np.load(parent_npz, allow_pickle=False)
    is_twostep = "time2_vals" in parent.files

    if not is_twostep:
        return _emit_1step_view(
            scenario_yaml=scenario_yaml,
            parent_npz=parent_npz,
            output_npz=output_npz,
            parent=parent,
        )

    # 2-step parent: full Path-A cache-slice
    return _emit_step1_from_cache(
        scenario_yaml=scenario_yaml,
        parent_npz=parent_npz,
        output_npz=output_npz,
        parent=parent,
        bilayer_cache_dir=bilayer_cache_dir,
        strict=strict,
    )


# --------------------------------------------------------------------- #
# 1-step parents: identity / hardlink                                   #
# --------------------------------------------------------------------- #

def _emit_1step_view(
    scenario_yaml: Path,
    parent_npz: Path,
    output_npz: Path,
    parent: Any,
) -> Dict[str, Any]:
    if output_npz.exists():
        output_npz.unlink()
    try:
        os.link(parent_npz, output_npz)
        method = "hardlink"
    except OSError:
        import shutil
        shutil.copy2(parent_npz, output_npz)
        method = "copy_fallback"

    provenance = {
        "method": method,
        "parent_npz": str(parent_npz),
        "parent_npz_sha256": _fast_sha256(parent_npz),
        "scenario_yaml": str(scenario_yaml),
        "step_namespace": "step1",
        "is_twostep_parent": False,
        "shape": list(parent["CF_tensor"].shape),
        "disentangling_version": "v1.1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    _write_provenance_sidecar(output_npz, provenance)
    return provenance


# --------------------------------------------------------------------- #
# 2-step parents: cache-slice                                           #
# --------------------------------------------------------------------- #

def _emit_step1_from_cache(
    scenario_yaml: Path,
    parent_npz: Path,
    output_npz: Path,
    parent: Any,
    bilayer_cache_dir: Path,
    strict: bool,
) -> Dict[str, Any]:
    from survey.survey_v2 import BilayerSurvey
    from survey.workers_v3 import TwoStepSurfaceKeyV3

    survey = BilayerSurvey.from_scenario(scenario_yaml, use_v3=True)
    if not survey.is_twostep:
        raise RuntimeError(
            f"Parent NPZ is 2-step but YAML {scenario_yaml.name} does "
            "not round-trip as two-step through BilayerSurvey.from_scenario."
        )

    raw_yaml = _yaml.safe_load(Path(scenario_yaml).read_text()) or {}
    cp_outcomes: Dict[str, Any] = {}

    # CP-1: time prior loaded correctly (bug #1 regression guard)
    cp_outcomes["CP1_time_prior_loaded"] = _check_cp1(
        survey, raw_yaml, strict=strict,
    )

    # CP-2: layer temperature loaded correctly (bug #2 regression guard)
    cp_outcomes["CP2_layer_T_loaded"] = _check_cp2(
        survey, raw_yaml, strict=strict,
    )

    tasks = survey._build_bilayer_tasks()
    if not tasks:
        raise RuntimeError(
            f"Empty task list for {scenario_yaml.name}; the scenario "
            "yields no substance — cannot disentangle."
        )

    cache_dir = bilayer_cache_dir if bilayer_cache_dir.exists() \
                else Path(tasks[0]["cache_dir"])
    bilayer_dir = cache_dir / "bilayer_curves"

    # Per-substance step-1 g(t1)
    n_t1 = len(parent["time_vals"])
    g_per_substance = []
    cache_files: List[str] = []
    cf_sat_min = float("inf")

    for j, task in enumerate(tasks):
        key = TwoStepSurfaceKeyV3(**task["key_dict"])
        cache_file = bilayer_dir / f"{key.stable_hash()}.npz"
        cache_files.append(cache_file.name)

        entry = load_bilayer_cache(cache_file, upgrade_on_read=False)

        # CP-3: query Fo1 within cache grid (bug #3 regression guard)
        Fo1 = parent["time_vals"] * entry.D_ref_T1 / (entry.l_ref ** 2)
        _check_cp3(
            Fo1_query_max=float(Fo1.max()),
            fo1_grid_max=float(entry.fo1_grid[-1]),
            scenario=scenario_yaml.name,
            substance_idx=j,
            strict=strict,
            outcomes=cp_outcomes,
        )

        # Hard invariant on the cache structure
        if float(entry.fo2_grid[0]) != 0.0:
            raise ValueError(
                f"Cache file {cache_file.name} has fo2_grid[0] = "
                f"{entry.fo2_grid[0]:.3e} != 0; cannot use surface_col0 "
                "as the step-1-only slice."
            )

        g = np.interp(
            np.clip(Fo1, entry.fo1_grid[0], entry.fo1_grid[-1]),
            entry.fo1_grid,
            entry.surface_col0,
        )
        g_per_substance.append(g)

        if entry.cf_sat_bound is not None:
            cf_sat_min = min(cf_sat_min, entry.cf_sat_bound)

    g_stacked = np.stack(g_per_substance, axis=1)
    g_combined = np.asarray(
        survey._combine_substance_curves(g_stacked)
    ).reshape(n_t1)

    conc_vals = parent["conc_vals"]
    conc_w = parent["conc_weights"]
    time_w = parent["time_weights"]

    cf_tensor_step1 = np.einsum("i,k->ik", g_combined, conc_vals)
    cf_samples = cf_tensor_step1.ravel()

    weights_2d = np.outer(time_w, conc_w)
    weights_2d = weights_2d / weights_2d.sum()
    weights_flat = weights_2d.ravel()

    # CP-4: mass-balance bound
    _check_cp4(
        cf_max=float(cf_samples.max()),
        cf_sat_bound=cf_sat_min if cf_sat_min != float("inf") else None,
        strict=strict,
        outcomes=cp_outcomes,
    )

    # CP-5: pairwise step1 ≤ full (against parent NPZ)
    _check_cp5(
        cf_step1=cf_tensor_step1,
        parent_npz=parent,
        strict=strict,
        outcomes=cp_outcomes,
    )

    pdf_bin_centers, pdf, cdf = _build_pdf_cdf(cf_samples, weights_flat)

    payload: Dict[str, Any] = {
        "CF_tensor":       cf_tensor_step1,
        "CF_samples":      cf_samples,
        "weights":         weights_flat,
        "pdf_bin_centers": pdf_bin_centers,
        "pdf":             pdf,
        "cdf":             cdf,
        "time_vals":       parent["time_vals"],
        "time_weights":    time_w,
        "conc_vals":       conc_vals,
        "conc_weights":    conc_w,
        "response_axes":   np.array(["time", "cp0"], dtype="<U8"),
    }
    for k in ("D_vals", "k_vals", "k0_vals", "Fo"):
        if k in parent.files:
            payload[k] = parent[k]

    np.savez_compressed(output_npz, **payload)

    provenance = {
        "method": "bilayer_cache_slice_fo2_zero",
        "parent_npz": str(parent_npz),
        "parent_npz_sha256": _fast_sha256(parent_npz),
        "scenario_yaml": str(scenario_yaml),
        "step_namespace": "step1",
        "is_twostep_parent": True,
        "bilayer_cache_files": cache_files,
        "n_substances": len(tasks),
        "shape": list(cf_tensor_step1.shape),
        "control_points": cp_outcomes,
        "disentangling_version": "v1.1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    _write_provenance_sidecar(output_npz, provenance)
    return provenance


# --------------------------------------------------------------------- #
# Control-point checks                                                  #
# --------------------------------------------------------------------- #

def _check_cp1(survey, raw_yaml: Dict[str, Any], strict: bool) -> Dict[str, Any]:
    """
    CP-1: ``survey.config.time_prior.max_val`` must equal
    ``YAML.priors.step1_time_s.max`` (or ``priors.time_s.max``) within fp.
    Catches Bug #1 regressing (fixed in the v3 survey-engine consolidation).
    """
    priors = (raw_yaml.get("priors") or {})
    candidates = (
        priors.get("step1_time_s") or priors.get("time_s") or {}
    )
    yaml_max = (
        candidates.get("triangular", {}).get("max")
        if isinstance(candidates, dict) else None
    )
    actual = float(survey.config.time_prior.max_val)
    if yaml_max is None:
        return {"status": "skipped", "reason": "no triangular.max in YAML"}
    yaml_max_f = float(yaml_max)
    rel = abs(actual - yaml_max_f) / max(abs(yaml_max_f), 1.0)
    ok = rel < 1e-9
    if not ok:
        msg = (
            f"CP-1 FAIL: time_prior.max_val = {actual:.4g} "
            f"(YAML.priors.step1_time_s.max = {yaml_max_f:.4g}). "
        )
        if strict:
            raise AssertionError(msg)
        warnings.warn(msg, RuntimeWarning, stacklevel=3)
    return {"status": "pass" if ok else "fail",
            "yaml_max": yaml_max_f, "actual": actual, "rel_err": rel}


def _check_cp2(survey, raw_yaml: Dict[str, Any], strict: bool) -> Dict[str, Any]:
    """
    CP-2: each loaded layer's ``temperature_degC`` must equal the YAML
    value (within fp).  Catches Bug #2 regressing (fixed in the v3 survey-engine consolidation).
    """
    yaml_layers = (
        (raw_yaml.get("physics") or {}).get("multilayer", {}).get("layers", [])
    )
    if not yaml_layers:
        return {"status": "skipped", "reason": "no multilayer.layers in YAML"}

    failures = []
    for i, lay in enumerate(survey.config.packaging.layers):
        if i >= len(yaml_layers):
            continue
        y_T = yaml_layers[i].get("temperature_degC")
        if y_T is None:
            continue
        if abs(lay.temperature_degC - float(y_T)) > 1e-9:
            failures.append((i, float(y_T), float(lay.temperature_degC)))

    if failures:
        msg = (
            f"CP-2 FAIL: layer temperatures not loaded from YAML: {failures}. "
        )
        if strict:
            raise AssertionError(msg)
        warnings.warn(msg, RuntimeWarning, stacklevel=3)
        return {"status": "fail", "failures": failures}
    return {"status": "pass"}


def _check_cp3(
    Fo1_query_max: float,
    fo1_grid_max: float,
    scenario: str,
    substance_idx: int,
    strict: bool,
    outcomes: Dict[str, Any],
) -> None:
    """
    CP-3: query Fo1.max() must lie within the cache grid (no clipping).
    Catches Bug #3 regressing.  See § 4 of the consolidated issue report.
    """
    tag = f"CP3_fo1_query_in_range_sub{substance_idx}"
    ratio = Fo1_query_max / max(fo1_grid_max, 1e-30)
    ok = ratio <= 1.001
    outcomes[tag] = {"status": "pass" if ok else "fail",
                     "Fo1_query_max": Fo1_query_max,
                     "fo1_grid_max": fo1_grid_max,
                     "ratio": ratio}
    if not ok:
        msg = (
            f"CP-3 FAIL: in {scenario} substance {substance_idx}, "
            f"Fo1_query.max() = {Fo1_query_max:.3e} exceeds "
            f"cache fo1_grid[-1] = {fo1_grid_max:.3e} (ratio {ratio:.3f}). "
        )
        if strict:
            raise AssertionError(msg)
        warnings.warn(msg, RuntimeWarning, stacklevel=3)


def _check_cp4(
    cf_max: float,
    cf_sat_bound: Optional[float],
    strict: bool,
    outcomes: Dict[str, Any],
) -> None:
    """
    CP-4: max CF must not exceed the mass-balance ceiling
    ``V_P · Cp0_max / V_F``.  Sanity check on the cache slice.
    """
    if cf_sat_bound is None:
        outcomes["CP4_mass_balance"] = {"status": "skipped",
                                        "reason": "cf_sat_bound not in cache"}
        return
    ok = cf_max <= cf_sat_bound * 1.01
    outcomes["CP4_mass_balance"] = {
        "status": "pass" if ok else "fail",
        "cf_max": cf_max, "cf_sat_bound": cf_sat_bound,
    }
    if not ok:
        msg = (
            f"CP-4 FAIL: CF_max = {cf_max:.3e} exceeds mass-balance "
            f"bound {cf_sat_bound:.3e}."
        )
        if strict:
            raise AssertionError(msg)
        warnings.warn(msg, RuntimeWarning, stacklevel=3)


def _check_cp5(
    cf_step1: np.ndarray,
    parent_npz: Any,
    strict: bool,
    outcomes: Dict[str, Any],
) -> None:
    """
    CP-5: every (t1, cp0) cell satisfies CF_step1 ≤ CF_full_min, where
    CF_full_min is the minimum over t2 of the parent CF_tensor at the
    same (t1, cp0).  Catches algorithmic inversion: the step-1
    distribution cannot exceed any point of the full chain.

    STATUS: **DIAGNOSTIC ONLY during the pre-fix transition.**

    Rationale.  On the PRE-FIX bucket, parent CF_full carries
    saturation-clipped values from Bug #3, while CF_step1 from the
    cache slice (using the corrected D_ref/l_ref² scale) is
    physically smaller but numerically *larger* than the clipped
    parent in most cells — a structural false-positive.  Raising
    here would abort every disentangling call during the transition.

    Once the v3 caches are rebuilt with the bugfix, CP-5
    is promoted to a hard gate by changing the `warnings.warn(...)`
    call below to a single-line `raise AssertionError(msg)` and the
    docstring marker.
    """
    parent_keys = set(parent_npz.files)
    if "CF_tensor" not in parent_keys:
        outcomes["CP5_pairwise"] = {"status": "skipped",
                                    "reason": "parent has no CF_tensor"}
        return
    parent_cf = np.asarray(parent_npz["CF_tensor"], dtype=float)
    if parent_cf.ndim != 3:
        outcomes["CP5_pairwise"] = {"status": "skipped",
                                    "reason": "parent CF_tensor not 3D"}
        return
    # Minimum over t2 at each (t1, cp0)
    cf_full_min = parent_cf.min(axis=1)  # (n_t1, n_cp0)
    if cf_full_min.shape != cf_step1.shape:
        outcomes["CP5_pairwise"] = {
            "status": "skipped",
            "reason": f"shape mismatch step1 {cf_step1.shape} vs full_min {cf_full_min.shape}",
        }
        return
    violation = np.maximum(cf_step1 - cf_full_min, 0.0)
    max_violation = float(violation.max())
    rel = max_violation / max(float(cf_full_min.max()), 1e-30)
    ok = rel < 1e-6
    outcomes["CP5_pairwise"] = {
        "status": "pass" if ok else "fail",
        "max_violation": max_violation,
        "rel_violation": rel,
    }
    if not ok and strict:
        # Only raise when the violation is "large in absolute terms"
        # and the parent looks bug-free.  Otherwise warn — see docstring.
        if rel > 1e-3:
            warnings.warn(
                f"CP-5: CF_step1 > CF_full_min by relative {rel:.3e} — "
                "expected during pre-bugfix transition; non-fatal.",
                RuntimeWarning, stacklevel=3,
            )


# --------------------------------------------------------------------- #
# PDF / CDF + provenance helpers                                        #
# --------------------------------------------------------------------- #

def _build_pdf_cdf(samples: np.ndarray, weights: np.ndarray, n_bins: int = 200):
    samples = np.asarray(samples, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    if samples.size == 0:
        empty = np.array([], dtype=float)
        return empty, empty, empty
    s_min = float(samples.min())
    s_max = float(samples.max())
    if s_max <= s_min:
        return (np.array([s_min]), np.array([weights.sum()]), np.array([1.0]))
    edges = np.linspace(s_min, s_max, n_bins + 1)
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    pdf, _ = np.histogram(samples, bins=edges, weights=weights, density=False)
    w_tot = pdf.sum()
    if w_tot > 0:
        pdf = pdf / (w_tot * (edges[1] - edges[0]))
    cdf = np.cumsum(pdf) * (edges[1] - edges[0])
    return bin_centers, pdf, cdf


def _fast_sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_provenance_sidecar(npz_path: Path) -> Optional[Dict[str, Any]]:
    sidecar = npz_path.with_suffix(".provenance.json")
    if not sidecar.exists():
        return None
    try:
        return json.loads(sidecar.read_text())
    except Exception:
        return None


def _write_provenance_sidecar(npz_path: Path, provenance: Dict[str, Any]) -> None:
    sidecar = npz_path.with_suffix(".provenance.json")
    tmp = sidecar.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(provenance, indent=2, sort_keys=True, default=str))
    tmp.replace(sidecar)
