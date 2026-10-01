"""
survey/workers_v3.py — optimised bilayer / two-step workers.

Relative to workers_v2:

- Two-step solver: per-t1 chained senspatankar + a single multi-t2 resume
  call. nt1 × nt2 inner solves collapse to 2 × nt1 solves (one step-1
  solve + one multi-t2 step-2 resume per t1 node).
- Profile propagation via sol1.resume(..., Cx0=previousCx): removes the
  "C0 = C0_initial" approximation in workers_v2 (line 294) — step 2 now
  starts from the real end-of-step-1 polymer profile.
- Equilibrium short-circuit on both axes using sol.PR.peq / sol.PR.k0
  (closed-form CFeq, computed by the solver from input geometry; free).

One-step solver: already optimal in workers_v2 (one senspatankar call
with t = fo_grid, multi-time output). workers_v3 wraps it unchanged but
clamps any trailing points with `|CF - CFeq| < rel_tol` to CFeq exactly
to remove late-time numerical drift near saturation.

Cache keys are versioned via the `worker_version` field so v2 and v3
cache entries never collide on disk.

@project: SFPPy — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np

from survey.workers import (
    make_fo_grid,
    assert_bilayer_mass_conservation,
    _cfeq_from_sol,
    _eq_knee,
)


# =====================================================================
# CFeq helper — canonical definitions live in survey.workers
# (_cfeq_from_sol, _eq_knee are re-exported above for backward
# compatibility with code that imports them from this module).
# =====================================================================


# =====================================================================
# Cache keys — versioned to avoid v2 collisions
# =====================================================================

@dataclass(frozen=True)
class BilayerCurveKeyV3:
    polymer_1: str
    D_1: float
    k_1: float
    l_1_m: float
    C0_1: float
    polymer_2: str
    D_2: float
    k_2: float
    l_2_m: float
    C0_2: float
    k0: float
    h: float
    surface_area: float
    food_volume: float
    contact_temperature_degC: float
    CF0: float
    Fo_max: float
    n_fo: int
    focal_layer: int
    # 'v3.1': food-volume normalisation fix in _normalise_bilayer
    # (V_ref = l_ref·A, was (l_1+l_2)·A). Bumped to invalidate the
    # mass-non-conserving curves cached under 'v3'.
    worker_version: str = 'v3.1'

    def stable_hash(self) -> str:
        raw = json.dumps(asdict(self), sort_keys=True, default=str)
        return hashlib.sha256(raw.encode()).hexdigest()[:16]


@dataclass(frozen=True)
class TwoStepSurfaceKeyV3:
    polymer_1: str
    D_1_T1: float
    D_1_T2: float
    k_1: float
    l_1_m: float
    C0_1: float
    polymer_2: str
    D_2_T1: float
    D_2_T2: float
    k_2: float
    l_2_m: float
    C0_2: float
    k0: float
    h: float
    surface_area: float
    food_volume: float
    T1_degC: float
    T2_degC: float
    CF0: float
    Fo1_max: float
    Fo2_max: float
    n_fo1: int
    n_fo2: int
    focal_layer: int
    # 'v3.1': food-volume normalisation fix in _normalise_bilayer
    # (V_ref = l_ref·A, was (l_1+l_2)·A).
    # 'v3.2': real-units bilayer chain (fix of the workers_v3.py:372 step-2
    # grid over-scaling). Bumped
    # to invalidate the over-integrated surfaces cached under 'v3.1'.
    worker_version: str = 'v3.2'

    def stable_hash(self) -> str:
        raw = json.dumps(asdict(self), sort_keys=True, default=str)
        return hashlib.sha256(raw.encode()).hexdigest()[:16]


@dataclass(frozen=True)
class MonolayerTwoStepKeyV3:
    """Cache key for a NATIVE monolayer two-step master surface g(Fo1,Fo2).

    A monolayer two-step is one layer, two temperatures, two chained time
    priors — NOT a degenerate bilayer. Both Fourier axes use a SINGLE common
    clock anchored to ``D_ref = max(D_T1, D_T2)`` (``reference_policy``), so the
    normalised diffusivity is ``<= 1`` in both steps. The slow (storage) step
    then has a large *formal* clock value but tiny ``D_norm = D_T1/D_ref`` —
    slow, non-stiff dynamics, NOT an enormous horizon with fast normalised
    dynamics; the fast (oven) step has ``D_norm = 1`` and no artificial
    cross-step ``×(D_T2/D_T1)`` horizon inflation.

    A distinct dataclass + ``worker_version`` gives a distinct ``stable_hash``
    space, so these entries can never collide with bilayer surfaces
    (``TwoStepSurfaceKeyV3``) and the shipped bilayer behaviour is untouched.
    Stored in the same ``bilayer_curves/`` cache dir and surface layout, so
    ``step_resolved_cf_tensors``, ``fo2=0`` slicing and the oven-effect
    analysis apply unchanged.
    """
    polymer: str
    D_T1: float
    D_T2: float
    k: float
    l_m: float
    C0: float
    k0: float
    h: float
    surface_area: float
    food_volume: float
    T1_degC: float
    T2_degC: float
    CF0: float
    Fo1_max: float
    Fo2_max: float
    n_fo1: int
    n_fo2: int
    # Explicit, auditable + cache-distinguishing reference policy.
    reference_policy: str = 'maxD_common_clock'
    # 'v3.2-mono': step-2 resume carries the storage-end food concentration
    #   (CF0=sol1.CFtarget) so CF(1+2) >= CF(1) and the fo2=0 slice == CF(1).
    # 'v3.3-mono': √Fo-spaced grids (uniform in √Fo) + higher resolution +
    #   timescale='sqrt' solver stepping. √-spacing concentrates points at high
    #   Fo on the max-D-stretched storage axis, giving ~1% interpolation
    #   accuracy where the t-prior actually lives. (Bilayer stays linear.)
    worker_version: str = 'v3.3-mono'

    def stable_hash(self) -> str:
        raw = json.dumps(asdict(self), sort_keys=True, default=str)
        return hashlib.sha256(raw.encode()).hexdigest()[:16]


# =====================================================================
# Normalisation helper (shared with v2 semantics)
# =====================================================================

def _normalise_bilayer(
    *, l_1_m, D_1, l_2_m, D_2, k_1, k_2, h, food_volume, surface_area,
):
    """
    Return (l_1_norm, l_2_norm, D_1_norm, D_2_norm, l_ref, D_ref,
            food_volume_norm, h_norm).

    Reference layer is the one with highest resistance D/(k·l).
    Identical to workers_v2's inline normalisation so surfaces are
    comparable.
    """
    perm_1 = D_1 / (k_1 * l_1_m) if k_1 * l_1_m > 0 else float('inf')
    perm_2 = D_2 / (k_2 * l_2_m) if k_2 * l_2_m > 0 else float('inf')
    if perm_1 <= perm_2:
        l_ref, D_ref = l_1_m, D_1
    else:
        l_ref, D_ref = l_2_m, D_2
    if D_ref <= 0:
        D_ref = max(D_1, D_2, 1e-30)

    l_1_norm = l_1_m / l_ref
    l_2_norm = l_2_m / l_ref
    D_1_norm = D_1 / D_ref
    D_2_norm = D_2 / D_ref

    # Food volume must be normalised on the SAME length scale as the layers
    # (l_ref·A), not the total packaging volume. The kernel's dilution number
    # is L = A·l_ref/V_F (migration.py:4168); with lengths divided by l_ref and
    # surfacearea=1, the food volume has to be V_F/(l_ref·A) for L_kernel to
    # equal L_real. Using (l_1+l_2)·A inflated L — and thus CF/CP0 — by exactly
    # (l_1+l_2)/l_ref, breaching mass conservation for bilayers (monolayers are
    # exempt since l_ref == l_total), as shown by a mass-bound falsification test.
    V_ref = l_ref * surface_area
    food_volume_norm = food_volume / V_ref if V_ref > 0 else food_volume
    h_norm = h * l_ref / D_ref if D_ref > 0 else h

    return (l_1_norm, l_2_norm, D_1_norm, D_2_norm,
            l_ref, D_ref, food_volume_norm, h_norm)


# =====================================================================
# 1-step bilayer master curve — v3 (equilibrium tail clamp)
# =====================================================================

def solve_bilayer_master_curve_v3(
    *,
    polymer_1: str, D_1: float, k_1: float, l_1_m: float, C0_1: float,
    polymer_2: str, D_2: float, k_2: float, l_2_m: float, C0_2: float,
    k0: float, h: float, surface_area: float, food_volume: float,
    contact_temperature_degC: float, CF0: float,
    Fo_max: float, n_fo: int,
    focal_layer: int = 0,
    fo_min_floor: float = 1e-15,
    eq_rel_tol: float = 1e-3,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    1-step master curve identical to workers_v2 with an equilibrium
    tail clamp. One senspatankar call (multi-t output); trailing
    samples within `eq_rel_tol` of CFeq are set to CFeq exactly.

    Returns (fo_grid, cf_over_cp0).
    """
    from patankar.layer import layer
    from patankar.food import foodphysics
    from patankar.migration import senspatankar

    (l_1_norm, l_2_norm, D_1_norm, D_2_norm,
     l_ref, D_ref, food_volume_norm, h_norm) = _normalise_bilayer(
        l_1_m=l_1_m, D_1=D_1, l_2_m=l_2_m, D_2=D_2,
        k_1=k_1, k_2=k_2, h=h,
        food_volume=food_volume, surface_area=surface_area,
    )

    T = contact_temperature_degC
    fo_grid = make_fo_grid(Fo_max=Fo_max, n_fo=n_fo, fo_min_floor=fo_min_floor)

    lay1 = layer(l=l_1_norm, D=D_1_norm, k=k_1, C0=C0_1, T=T)
    lay2 = layer(l=l_2_norm, D=D_2_norm, k=k_2, C0=C0_2, T=T)
    bilayer = lay1 + lay2

    food = foodphysics(
        k=k_1 if focal_layer == 0 else k_2,
        k0=k0, h=h_norm,
        surfacearea=1.0, volume=food_volume_norm,
        contacttime=float(fo_grid[-1]),
        contacttemperature=T, CF0=CF0,
    )

    sol = senspatankar(bilayer, food, t=fo_grid, autotime=False)
    # sol.t may be extended beyond fo_grid for post-contact diagnostics;
    # route through the shared helper (workplan § 1.6).
    from survey.utils import cf_at_user_grid
    cf = cf_at_user_grid(sol, fo_grid)

    cfeq = _cfeq_from_sol(sol)
    knee = _eq_knee(cf, cfeq, eq_rel_tol)
    if knee < len(cf):
        cf[knee:] = cfeq

    assert_bilayer_mass_conservation(
        cf, l_1_m=l_1_m, l_2_m=l_2_m, C0_1=C0_1, C0_2=C0_2,
        surface_area=surface_area, food_volume=food_volume, context="(1-step v3)",
    )
    return fo_grid, cf


# =====================================================================
# 2-step master surface — v3 (per-t1 chained, multi-t2 resume)
# =====================================================================

def solve_twostep_master_surface_v3(
    *,
    polymer_1: str, D_1_T1: float, D_1_T2: float,
    k_1: float, l_1_m: float, C0_1: float,
    polymer_2: str, D_2_T1: float, D_2_T2: float,
    k_2: float, l_2_m: float, C0_2: float,
    k0: float, h: float, surface_area: float, food_volume: float,
    T1_degC: float, T2_degC: float, CF0: float,
    Fo1_max: float, Fo2_max: float,
    n_fo1: int = 30, n_fo2: int = 15,
    focal_layer: int = 0,
    fo_min_floor: float = 1e-15,
    eq_rel_tol: float = 1e-3,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    2-step master surface — REAL-UNITS bilayer chain (fix 2026-06-21).

    Computed as the 2-layer analogue of ``mono_twostep_chained.chain_substance``:
    PHYSICAL diffusivities are passed to the engine and its INTERNAL Fourier
    normalisation does the scaling. The engine classifies long/short-term by the
    MAGNITUDE of the D it receives (the 1e-12/1e-6 m²/s thresholds in migration.py),
    so it must be given real D, never D/D_ref. For each fo1_i:
      (a) solve step 1 at T1 in real units up to the physical time
          ``t1_i = Fo1_i·l_ref²/D_ref_T1``;
      (b) resume step 2 at T2 (real D@T2) with the end-of-step-1 polymer profile,
          re-passing ``CF0 = sol1.CFtarget`` so the storage-end food mass is carried
          (else the fo2≈0 column reads ~0 and CF12<CF1), requesting CF at all fo2
          points in one solve.

    Axes are returned in the SAME dimensionless ``(fo1_grid, fo2_grid)`` coordinates
    the cache key and the query side use — only the SOLVE is in real units. ``fo1_grid``
    is in ``D_ref_T1`` units, ``fo2_grid`` in ``D_ref_T2`` units (matching ``Fo1_max``/
    ``Fo2_max`` from the key and the query scales ``α_T1=D_ref_T1/l_ref²``,
    ``α_T2=D_ref_T2/l_ref²``), so the t↔Fo round-trip is exact.

    Supersedes the normalised-surface implementation (``worker_version`` v3.1), whose
    step-2 solver grid was over-scaled by ``D_ref_T2/D_ref_T1`` under a false premise
    (the grid is already in T2-ref units) — that over-integrated the oven step toward
    equilibrium (~33× over-estimate on frozen→reheat) AND timed out.

    Equilibrium short-circuit on both axes:
      - Step 1 past CFeq_T1: freeze the last step-1 solution and reuse.
      - Step 2 past CFeq_T2 within a row: clamp trailing samples to CFeq_T2.

    Returns (fo1_grid, fo2_grid, surface).
    """
    # This is now a THIN ADAPTER over the unified real-units core
    # `survey.chain_engine.solve_chain_master_surface` (N-layer × K-step). The bilayer
    # two-step is the N=2, K=2 case. Verified BIT-IDENTICAL to the previous inline
    # implementation (T8_engine_regression: max_rel/abs_diff = 0.0). The bilayer-specific
    # mass-conservation assert is preserved here.
    from survey.chain_engine import LayerSpec, solve_chain_master_surface

    layers = [
        LayerSpec(l_m=l_1_m, D_T1=D_1_T1, k=k_1, C0=C0_1, D_T2=D_1_T2),
        LayerSpec(l_m=l_2_m, D_T1=D_2_T1, k=k_2, C0=C0_2, D_T2=D_2_T2),
    ]
    fo1_grid, fo2_grid, surface = solve_chain_master_surface(
        layers=layers, k0=k0, h=h, surface_area=surface_area, food_volume=food_volume,
        T1_degC=T1_degC, CF0=CF0, Fo1_max=Fo1_max, n_fo1=n_fo1, focal_layer=focal_layer,
        T2_degC=T2_degC, Fo2_max=Fo2_max, n_fo2=n_fo2,
        fo_min_floor=fo_min_floor, eq_rel_tol=eq_rel_tol,
    )

    assert_bilayer_mass_conservation(
        surface, l_1_m=l_1_m, l_2_m=l_2_m, C0_1=C0_1, C0_2=C0_2,
        surface_area=surface_area, food_volume=food_volume,
        context="(2-step v3 real-units, chain_engine adapter)",
    )
    return fo1_grid, fo2_grid, surface


# =====================================================================
# n>=3 spatial layers — generalisation
#
# The solve core (`survey.chain_engine.solve_chain_master_surface`) is already
# N-layer × K-step. What was missing is the survey-side arity: cache
# keys, workers and task building were hardcoded to 2 layers (six sites). The n<=2
# classes above are deliberately UNTOUCHED (byte-identical keys → the warm
# g-cache and all P3-validated bilayer paths are preserved); n>=3 gets its own
# tuple-valued keys in a distinct `worker_version` hash space.
# =====================================================================

@dataclass(frozen=True)
class MultilayerCurveKeyV3:
    """Cache key for an n>=3-layer SINGLE-STEP master curve (real-units chain).

    Per-layer parameters are tuples ordered food-contact first (index 0 =
    contact layer, matching the scenario layer order). Solved by
    `chain_engine.solve_chain_master_surface` (K=1) in REAL units — the
    normalised 1-step bilayer path (`_normalise_bilayer`) is not used for
    n>=3. Distinct dataclass + worker_version → distinct stable_hash space:
    can never collide with Bilayer/TwoStep/Monolayer keys.
    """
    polymers: Tuple[str, ...]
    D_T1s: Tuple[float, ...]
    ks: Tuple[float, ...]
    l_ms: Tuple[float, ...]
    C0s: Tuple[float, ...]
    k0: float
    h: float
    surface_area: float
    food_volume: float
    contact_temperature_degC: float
    CF0: float
    Fo_max: float
    n_fo: int
    focal_layer: int
    worker_version: str = 'v3.3-multi'

    def stable_hash(self) -> str:
        raw = json.dumps(asdict(self), sort_keys=True, default=str)
        return hashlib.sha256(raw.encode()).hexdigest()[:16]


@dataclass(frozen=True)
class MultilayerSurfaceKeyV3:
    """Cache key for an n>=3-layer TWO-STEP master surface g(Fo1,Fo2).

    Same tuple layout as MultilayerCurveKeyV3 plus per-layer D at T2. Real-units
    chain (per-t1 solve + multi-t2 resume) via `chain_engine` — the same proven
    algorithm as the bilayer v3.2 surface, layer count generalised.
    """
    polymers: Tuple[str, ...]
    D_T1s: Tuple[float, ...]
    D_T2s: Tuple[float, ...]
    ks: Tuple[float, ...]
    l_ms: Tuple[float, ...]
    C0s: Tuple[float, ...]
    k0: float
    h: float
    surface_area: float
    food_volume: float
    T1_degC: float
    T2_degC: float
    CF0: float
    Fo1_max: float
    Fo2_max: float
    n_fo1: int
    n_fo2: int
    focal_layer: int
    worker_version: str = 'v3.3-multi'

    def stable_hash(self) -> str:
        raw = json.dumps(asdict(self), sort_keys=True, default=str)
        return hashlib.sha256(raw.encode()).hexdigest()[:16]


def assert_multilayer_mass_conservation(
    field: np.ndarray, *, l_ms, C0s, surface_area: float, food_volume: float,
    context: str = "",
) -> None:
    """n-layer mass-balance guard: max(CF/CP0) <= sum_i(C0_i·l_i)·A / V_F.

    The n-layer analogue of `assert_bilayer_mass_conservation` (same physical
    ceiling, same tolerance) — non-negotiable at solve time.
    """
    from survey.workers import MASS_CONSERVATION_TOL
    m0 = sum(float(c) * float(l) for c, l in zip(C0s, l_ms)) * surface_area
    if food_volume <= 0 or m0 <= 0:
        return
    ceiling = m0 / food_volume
    mx = float(np.nanmax(field))
    if mx > ceiling * (1.0 + MASS_CONSERVATION_TOL):
        raise ValueError(
            f"multilayer mass bound violated {context}: max CF/CP0={mx:.6g} "
            f"> ceiling={ceiling:.6g} (m0/V_F, {len(l_ms)} layers)")


def solve_multilayer_master_curve_v3(
    *, polymers, D_T1s, ks, l_ms, C0s,
    k0: float, h: float, surface_area: float, food_volume: float,
    contact_temperature_degC: float, CF0: float,
    Fo_max: float, n_fo: int, focal_layer: int = 0,
    fo_min_floor: float = 1e-15, eq_rel_tol: float = 1e-3,
) -> Tuple[np.ndarray, np.ndarray]:
    """n>=3-layer 1-step master curve via the real-units chain engine (K=1)."""
    from survey.chain_engine import LayerSpec, solve_chain_master_surface
    layers = [LayerSpec(l_m=float(l), D_T1=float(D), k=float(k), C0=float(c))
              for l, D, k, c in zip(l_ms, D_T1s, ks, C0s)]
    fo_grid, cf = solve_chain_master_surface(
        layers=layers, k0=k0, h=h, surface_area=surface_area,
        food_volume=food_volume, T1_degC=contact_temperature_degC, CF0=CF0,
        Fo1_max=Fo_max, n_fo1=n_fo, focal_layer=focal_layer,
        fo_min_floor=fo_min_floor, eq_rel_tol=eq_rel_tol,
    )
    assert_multilayer_mass_conservation(
        cf, l_ms=l_ms, C0s=C0s, surface_area=surface_area,
        food_volume=food_volume, context="(1-step v3.3-multi)")
    return fo_grid, cf


def solve_multilayer_surface_v3(
    *, polymers, D_T1s, D_T2s, ks, l_ms, C0s,
    k0: float, h: float, surface_area: float, food_volume: float,
    T1_degC: float, T2_degC: float, CF0: float,
    Fo1_max: float, Fo2_max: float, n_fo1: int = 30, n_fo2: int = 15,
    focal_layer: int = 0,
    fo_min_floor: float = 1e-15, eq_rel_tol: float = 1e-3,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """n>=3-layer 2-step master surface via the real-units chain engine (K=2)."""
    from survey.chain_engine import LayerSpec, solve_chain_master_surface
    layers = [LayerSpec(l_m=float(l), D_T1=float(D1), k=float(k), C0=float(c),
                        D_T2=float(D2))
              for l, D1, D2, k, c in zip(l_ms, D_T1s, D_T2s, ks, C0s)]
    fo1_grid, fo2_grid, surface = solve_chain_master_surface(
        layers=layers, k0=k0, h=h, surface_area=surface_area,
        food_volume=food_volume, T1_degC=T1_degC, CF0=CF0,
        Fo1_max=Fo1_max, n_fo1=n_fo1, focal_layer=focal_layer,
        T2_degC=T2_degC, Fo2_max=Fo2_max, n_fo2=n_fo2,
        fo_min_floor=fo_min_floor, eq_rel_tol=eq_rel_tol,
    )
    assert_multilayer_mass_conservation(
        surface, l_ms=l_ms, C0s=C0s, surface_area=surface_area,
        food_volume=food_volume, context="(2-step v3.3-multi)")
    return fo1_grid, fo2_grid, surface


# =====================================================================
# Native monolayer two-step master surface (max-D common clock)
# =====================================================================

def _sqrt_fo_grid(Fo_max: float, n_fo: int) -> np.ndarray:
    """√Fo-uniform grid on ``[0, Fo_max]`` (includes 0).

    Uniform in √Fo, i.e. ``Fo = s²`` for ``s`` linearly spaced in
    ``[0, √Fo_max]``. Because a 1-D diffusion master curve ``CF/CP0`` is
    ≈ linear in √Fo, √-uniform sampling gives near-uniform interpolation
    accuracy and (on the max-D-stretched storage axis) places most points at
    high Fo where the t-prior lives — unlike log-spacing, which wastes points
    at tiny Fo.
    """
    Fo_max = float(max(Fo_max, 0.0))
    if Fo_max <= 0.0:
        return np.array([0.0, 1.0], dtype=float)
    s = np.linspace(0.0, np.sqrt(Fo_max), int(max(4, n_fo)))
    return np.unique(s * s)


def solve_monolayer_twostep_surface_v3(
    *,
    polymer: str, D_T1: float, D_T2: float,
    k: float, l_m: float, C0: float,
    k0: float, h: float, surface_area: float, food_volume: float,
    T1_degC: float, T2_degC: float, CF0: float,
    Fo1_max: float, Fo2_max: float,
    n_fo1: int = 30, n_fo2: int = 15,
    reference_policy: str = 'maxD_common_clock',
    fo_min_floor: float = 1e-15,
    eq_rel_tol: float = 1e-3,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Native monolayer two-step master surface ``g(Fo1,Fo2) = CF/CP0``.

    One layer chained ``T1 -> T2`` (storage then oven). A SINGLE common clock
    anchored to ``D_ref = max(D_T1, D_T2)``:

      - both Fo axes are ``D_ref·t/l²``;
      - each step layer carries normalised diffusivity ``D_step/D_ref <= 1``;
      - ``h_norm = h·l/D_ref`` preserves the physical Biot through
        ``h_norm/D_norm = h·l/D_step``.

    The slow (storage) step has a large *formal* clock horizon but tiny
    ``D_norm`` → slow, **non-stiff** dynamics (the numerical win is that the
    fast step is no longer an enormous horizon with fast normalised dynamics,
    as it was under the step-1-anchored ``×(D_T2/D_T1)`` rescale). The oven
    step has ``D_norm = 1`` and no artificial horizon inflation.

    Surface layout ``(fo1_grid, fo2_grid, surface)`` and equilibrium clamps are
    identical to the bilayer two-step, so the downstream slicing/analysis is
    unchanged. This is a *real one-layer* solve — not a degenerate bilayer.
    """
    from patankar.layer import layer
    from patankar.food import foodphysics
    from patankar.migration import senspatankar
    from survey.utils import cf_at_user_grid

    D_ref = max(float(D_T1), float(D_T2))
    if D_ref <= 0:
        D_ref = max(float(D_T1), float(D_T2), 1e-30)
    l_ref = l_m
    D_norm_T1 = D_T1 / D_ref          # <= 1
    D_norm_T2 = D_T2 / D_ref          # <= 1 (one of the two == 1)

    # Food volume on the l_ref·A scale (monolayer: l_ref == l). h on D_ref so
    # the physical Biot h·l/D_step is preserved via h_norm/D_norm.
    V_ref = l_ref * surface_area
    food_volume_norm = food_volume / V_ref if V_ref > 0 else food_volume
    h_norm = h * l_ref / D_ref if D_ref > 0 else h

    # √Fo-uniform grids (CF/CP0 is ≈ linear in √Fo, so √-uniform sampling
    # minimises linear-interpolation error and concentrates points at high Fo
    # on the max-D-stretched storage axis — unlike log-spacing). Includes 0.
    fo1_grid = _sqrt_fo_grid(Fo1_max, n_fo1)
    fo2_grid = _sqrt_fo_grid(Fo2_max, n_fo2)
    # Single common clock → step-2 solver time grid IS fo2_grid (no D-ratio
    # rescale, unlike the bilayer solver). senspatankar needs a strictly
    # increasing t>0 grid; the leading 0 is replaced by a genuinely tiny value
    # (fo_min_floor ≪ fo2_grid[1]) so column 0 samples fo2≈0 (storage alone),
    # preserving the fo2=0 slice semantics.
    fo2_solve = fo2_grid.copy()
    if fo2_solve[0] <= 0:
        fo2_solve = np.concatenate(([fo_min_floor], fo2_solve[1:]))

    lay_T1 = layer(l=1.0, D=D_norm_T1, k=k, C0=C0, T=T1_degC)
    lay_T2 = layer(l=1.0, D=D_norm_T2, k=k, C0=C0, T=T2_degC)

    food_T2 = foodphysics(
        k=k, k0=k0, h=h_norm, surfacearea=1.0, volume=food_volume_norm,
        contacttime=float(fo2_solve[-1]), contacttemperature=T2_degC, CF0=CF0,
    )

    surface = np.empty((len(fo1_grid), len(fo2_grid)), dtype=float)
    cfeq_T1 = None
    cfeq_T2 = None
    sol1_frozen = None

    for i, Fo1 in enumerate(fo1_grid):
        Fo1_use = float(Fo1) if Fo1 > 0 else fo_min_floor
        if sol1_frozen is not None:
            sol1 = sol1_frozen
        else:
            food_T1 = foodphysics(
                k=k, k0=k0, h=h_norm, surfacearea=1.0, volume=food_volume_norm,
                contacttime=Fo1_use, contacttemperature=T1_degC, CF0=CF0,
            )
            # timescale='sqrt' → solver steps uniformly in √Fo (diffusion-
            # natural); inherited by the step-2 resume via inputs["timescale"].
            sol1 = senspatankar(lay_T1, food_T1,
                                t=np.array([0.0, Fo1_use]), autotime=False,
                                timescale='sqrt')
            cf_end_step1 = float(sol1.CFtarget)
            if cfeq_T1 is None:
                cfeq_T1 = _cfeq_from_sol(sol1)
            if cfeq_T1 > 0 and abs(cf_end_step1 - cfeq_T1) <= eq_rel_tol * cfeq_T1:
                sol1_frozen = sol1

        # Chain into the oven step carrying BOTH the depleted polymer profile
        # (Cxprevious, inherited by resume) AND the storage-end food
        # concentration. resume() only auto-carries the food CF0 when the
        # medium is NOT overridden; since we override medium=food_T2 (T2, h@T2)
        # we must re-pass CF0=sol1.CFtarget, else the storage-migrated food mass
        # is discarded (fo2=0 slice would read ~0 instead of CF(1), and CF(1+2)
        # would fall below CF(1) — non-physical). See migration.resume:1953-1959.
        sol2 = sol1.resume(multilayer=lay_T2, medium=food_T2,
                           CF0=float(sol1.CFtarget),
                           t=fo2_solve, autotime=False)
        cf_row = cf_at_user_grid(sol2, fo2_solve)
        if cfeq_T2 is None:
            cfeq_T2 = _cfeq_from_sol(sol2)
        if cfeq_T2 is not None and cfeq_T2 > 0:
            knee = _eq_knee(cf_row, cfeq_T2, eq_rel_tol)
            if knee < len(cf_row):
                cf_row[knee:] = cfeq_T2
        surface[i, :] = cf_row

    # Mass guard: monolayer ceiling C0·l·A/V_F (l_2_m=0, C0_2=0 reduce the
    # shared helper to the one-layer source mass — not a bilayer representation).
    assert_bilayer_mass_conservation(
        surface, l_1_m=l_m, l_2_m=0.0, C0_1=C0, C0_2=0.0,
        surface_area=surface_area, food_volume=food_volume,
        context="(2-step monolayer v3, maxD clock)",
    )
    return fo1_grid, fo2_grid, surface


# =====================================================================
# Cache — reuses v2's BilayerCache directory layout
# =====================================================================

class BilayerCacheV3:
    """
    Persistent cache for v3 bilayer master curves / surfaces.

    Layout: cache_dir/bilayer_curves/<hash>.npz — same directory as v2;
    v3 vs v2 collisions are avoided by the `worker_version` field in the
    key (different hashes).
    """

    def __init__(self, cache_dir: str):
        self.root = Path(cache_dir).expanduser().resolve() / "bilayer_curves"
        self.root.mkdir(parents=True, exist_ok=True)
        self._stats = {"hits": 0, "misses": 0}

    def _path(self, key) -> Path:
        return self.root / f"{key.stable_hash()}.npz"

    def exists(self, key) -> bool:
        return self._path(key).exists()

    def load_1d(self, key: BilayerCurveKeyV3) -> Tuple[np.ndarray, np.ndarray]:
        data = np.load(self._path(key))
        self._stats["hits"] += 1
        return data["fo_grid"], data["cf_over_cp0"]

    def save_1d(self, key: BilayerCurveKeyV3, fo: np.ndarray, cf: np.ndarray):
        np.savez_compressed(self._path(key), fo_grid=fo, cf_over_cp0=cf)
        meta = self.root / f"{key.stable_hash()}.json"
        meta.write_text(json.dumps(asdict(key), indent=2, sort_keys=True, default=str))
        self._stats["misses"] += 1

    def load_2d(self, key: TwoStepSurfaceKeyV3) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        data = np.load(self._path(key))
        self._stats["hits"] += 1
        return data["fo1_grid"], data["fo2_grid"], data["surface"]

    def save_2d(self, key: TwoStepSurfaceKeyV3,
                fo1: np.ndarray, fo2: np.ndarray, surface: np.ndarray):
        np.savez_compressed(self._path(key), fo1_grid=fo1, fo2_grid=fo2, surface=surface)
        meta = self.root / f"{key.stable_hash()}.json"
        meta.write_text(json.dumps(asdict(key), indent=2, sort_keys=True, default=str))
        self._stats["misses"] += 1

    @property
    def stats(self) -> Dict[str, int]:
        return dict(self._stats)


# =====================================================================
# Worker entry points (parallel-safe)
# =====================================================================

def worker_bilayer_curve_v3(payload: Dict[str, Any]) -> Dict[str, Any]:
    key = BilayerCurveKeyV3(**payload["key_dict"])
    cache = BilayerCacheV3(payload["cache_dir"])
    if cache.exists(key):
        fo, cf = cache.load_1d(key)
        return {"status": "hit", "key_hash": key.stable_hash(), "fo": fo, "cf": cf}
    fo, cf = solve_bilayer_master_curve_v3(**{
        k: v for k, v in asdict(key).items()
        if k not in ('focal_layer', 'worker_version')
    }, focal_layer=key.focal_layer)
    cache.save_1d(key, fo, cf)
    return {"status": "miss", "key_hash": key.stable_hash(), "fo": fo, "cf": cf}


def worker_twostep_surface_v3(payload: Dict[str, Any]) -> Dict[str, Any]:
    key = TwoStepSurfaceKeyV3(**payload["key_dict"])
    cache = BilayerCacheV3(payload["cache_dir"])
    if cache.exists(key):
        fo1, fo2, surf = cache.load_2d(key)
        return {"status": "hit", "key_hash": key.stable_hash(),
                "fo1": fo1, "fo2": fo2, "surface": surf}
    fo1, fo2, surf = solve_twostep_master_surface_v3(**{
        k: v for k, v in asdict(key).items()
        if k not in ('focal_layer', 'worker_version')
    }, focal_layer=key.focal_layer)
    cache.save_2d(key, fo1, fo2, surf)
    return {"status": "miss", "key_hash": key.stable_hash(),
            "fo1": fo1, "fo2": fo2, "surface": surf}


def worker_multilayer_curve_v3(payload: Dict[str, Any]) -> Dict[str, Any]:
    """n>=3-layer 1-step curve worker. Same cache store, distinct
    MultilayerCurveKeyV3 hash space."""
    kd = dict(payload["key_dict"])
    for f in ("polymers", "D_T1s", "ks", "l_ms", "C0s"):
        kd[f] = tuple(kd[f])          # asdict/pickle round-trips may deliver lists
    key = MultilayerCurveKeyV3(**kd)
    cache = BilayerCacheV3(payload["cache_dir"])
    if cache.exists(key):
        fo, cf = cache.load_1d(key)
        return {"status": "hit", "key_hash": key.stable_hash(), "fo": fo, "cf": cf}
    fo, cf = solve_multilayer_master_curve_v3(**{
        k: v for k, v in asdict(key).items()
        if k not in ('focal_layer', 'worker_version')
    }, focal_layer=key.focal_layer)
    cache.save_1d(key, fo, cf)
    return {"status": "miss", "key_hash": key.stable_hash(), "fo": fo, "cf": cf}


def worker_multilayer_surface_v3(payload: Dict[str, Any]) -> Dict[str, Any]:
    """n>=3-layer 2-step surface worker. Same cache store, distinct
    MultilayerSurfaceKeyV3 hash space."""
    kd = dict(payload["key_dict"])
    for f in ("polymers", "D_T1s", "D_T2s", "ks", "l_ms", "C0s"):
        kd[f] = tuple(kd[f])
    key = MultilayerSurfaceKeyV3(**kd)
    cache = BilayerCacheV3(payload["cache_dir"])
    if cache.exists(key):
        fo1, fo2, surf = cache.load_2d(key)
        return {"status": "hit", "key_hash": key.stable_hash(),
                "fo1": fo1, "fo2": fo2, "surface": surf}
    fo1, fo2, surf = solve_multilayer_surface_v3(**{
        k: v for k, v in asdict(key).items()
        if k not in ('focal_layer', 'worker_version')
    }, focal_layer=key.focal_layer)
    cache.save_2d(key, fo1, fo2, surf)
    return {"status": "miss", "key_hash": key.stable_hash(),
            "fo1": fo1, "fo2": fo2, "surface": surf}


def worker_monolayer_twostep_surface_v3(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Native monolayer two-step surface worker (max-D common clock).

    Reuses the BilayerCacheV3 2-D store; the distinct MonolayerTwoStepKeyV3
    hash space keeps these entries separate from bilayer surfaces.
    """
    key = MonolayerTwoStepKeyV3(**payload["key_dict"])
    cache = BilayerCacheV3(payload["cache_dir"])
    if cache.exists(key):
        fo1, fo2, surf = cache.load_2d(key)
        return {"status": "hit", "key_hash": key.stable_hash(),
                "fo1": fo1, "fo2": fo2, "surface": surf}
    fo1, fo2, surf = solve_monolayer_twostep_surface_v3(**{
        k: v for k, v in asdict(key).items()
        if k not in ('worker_version',)
    })
    cache.save_2d(key, fo1, fo2, surf)
    return {"status": "miss", "key_hash": key.stable_hash(),
            "fo1": fo1, "fo2": fo2, "surface": surf}
