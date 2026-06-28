"""
survey/workers.py — Parallel Worker Functions for Master Curve Computation
==========================================================================

Worker functions compatible with multiprocessing for computing
master curves in parallel.

@project: SFPPy — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

import os
import math
import time
from typing import Dict, Any, Tuple

import numpy as np

# Keep OMP single-threaded to avoid oversubscription in workers
os.environ.setdefault("OMP_NUM_THREADS", "1")


# Default tolerance for the in-solver bilayer mass-conservation guard.
# The (l_1+l_2)/l_ref normalisation defect (now fixed) inflated CF/CP0 by
# >=1.18x; conserving solutions land at <=1.00x. 5% cleanly separates them.
MASS_CONSERVATION_TOL = 0.05


def assert_bilayer_mass_conservation(
    field: np.ndarray,
    *,
    l_1_m: float,
    l_2_m: float,
    C0_1: float,
    C0_2: float,
    surface_area: float,
    food_volume: float,
    tol: float = MASS_CONSERVATION_TOL,
    context: str = "",
) -> None:
    """In-solver guard: the dimensionless field ``max(CF/CP0)`` cannot exceed
    the source mass per unit food volume.

    ``CF/CP0`` is frame-independent (the kernel returns the physical ratio), so
    the guard MUST test it against the ceiling built from the REAL geometry —
    not the normalised one. The kernel always conserves mass in whatever frame
    it is handed; a wrong normalisation inflates BOTH the curve and a
    normalised ceiling in lock-step, so only the physical ceiling exposes it.

    Only the focal layer carries C0 (the other layer's C0 is 0), so the source
    mass is ``m0 = (C0_1*l_1_m + C0_2*l_2_m)*A`` and the absolute ceiling
    (no polymer retention) is ``m0 / V_F`` (migration.py:4186-4190).

    Raises ``ValueError`` on breach so a regression in the food-volume
    normalisation (``_normalise_bilayer``) fails loudly at solve time rather
    than silently emitting non-physical, mass-violating CF/CP0 curves.
    """
    m0 = (C0_1 * l_1_m + C0_2 * l_2_m) * surface_area
    if food_volume <= 0 or m0 <= 0:
        return
    ceiling = m0 / food_volume                       # V_source / V_F (physical)
    peak = float(np.nanmax(field))
    if peak > ceiling * (1.0 + tol):
        raise ValueError(
            f"bilayer mass-conservation breach"
            f"{(' ' + context) if context else ''}: "
            f"max(CF/CP0)={peak:.6g} > V_focal/V_F={ceiling:.6g} "
            f"(x{peak / ceiling:.3f}, tol={tol:.0%}). "
            f"Check _normalise_bilayer food-volume scaling (l_ref*A)."
        )


def make_fo_grid(Fo_max: float, n_fo: int, fo_min_floor: float = 1e-15) -> np.ndarray:
    """
    Create log-spaced Fourier number grid with explicit 0 at start.

    Parameters
    ----------
    Fo_max : float
        Maximum Fourier number.
    n_fo : int
        Number of grid points.
    fo_min_floor : float
        Floor value for log spacing.

    Returns
    -------
    np.ndarray
        Fo grid including 0.
    """
    Fo_max = float(max(Fo_max, 0.0))
    if Fo_max <= 0.0:
        return np.array([0.0, 1.0], dtype=float)

    # Fomin = max(1e-8, 1e-6·Fomax): grid floor is a fixed fraction of the
    # horizon (self-similar in Fo), not an absolute 1e-15 that over-resolves
    # the early region for large Fo. Mirrors the senspatankar first_step rule.
    fo_min_floor = max(fo_min_floor, 1e-8, 1e-6 * Fo_max)
    lo = math.log10(max(fo_min_floor, 1e-30))
    hi = math.log10(max(Fo_max, fo_min_floor))
    grid = np.logspace(lo, hi, int(max(4, n_fo)), dtype=float)
    grid = np.unique(np.concatenate([np.array([0.0], dtype=float), grid]))
    return grid


def _cfeq_from_sol(sol) -> float:
    """
    Equilibrium CF computed by the solver from input geometry.

    Line 4146-4149 of patankar/migration.py:
        peq  = k0 * m0P / (VF + sum((k0/k)*V))
        mFeq = (peq/k0) * VF

    CFeq = mFeq / VF is the mass-balance equilibrium; divided by VF gives
    concentration in the food. sol.PR exposes mFeq and VF directly.
    mFeq is an array (one entry per layer under multi-source initial
    conditions); the food-side equilibrium CF is the sum over layers
    divided by VF (total mass in food / VF).
    """
    mFeq = np.asarray(sol.PR.mFeq).reshape(-1)
    VF = float(np.asarray(sol.PR.VF).reshape(-1)[0])
    if VF <= 0:
        return 0.0
    return float(mFeq.sum() / VF)


def _eq_knee(cf_array: np.ndarray, cfeq: float, rel_tol: float) -> int:
    """
    Index of the first sample within rel_tol of equilibrium (inclusive).
    Returns len(cf_array) if equilibrium is never reached.
    """
    if cfeq <= 0:
        return len(cf_array)
    mask = np.abs(cf_array - cfeq) <= rel_tol * cfeq
    if not mask.any():
        return len(cf_array)
    return int(np.argmax(mask))


def solve_master_curve(
    *,
    k: float,
    k0: float,
    lP_m: float,
    D_real: float,
    Fo_max: float,
    n_fo: int,
    contact_temperature_degC: float,
    h: float,
    surface_area: float,
    food_volume: float,
    CF0: float,
    fo_min_floor: float = 1e-15,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Solve master curve g(Fo) = CF/CP0 using Patankar solver.

    Uses normalized coordinates (l=1, D=1, CP0=1) with time axis = Fo.
    Preserves physics through:
    1. Volume ratio: V_F_norm = V_F_real / V_P_real
    2. Biot number: h_norm = h_real * l_real / D_real

    Parameters
    ----------
    k : float
        Partition coefficient.
    k0 : float
        k0 parameter.
    lP_m : float
        Layer thickness (m).
    D_real : float
        Real diffusivity (m²/s).
    Fo_max : float
        Maximum Fourier number.
    n_fo : int
        Number of Fo grid points.
    contact_temperature_degC : float
        Contact temperature (°C).
    h : float
        Mass transfer coefficient (m/s).
    surface_area : float
        Contact surface area (m²).
    food_volume : float
        Food volume (m³).
    CF0 : float
        Initial food concentration.
    fo_min_floor : float
        Floor for Fo log grid.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (fo_grid, cf_over_cp0)
    """
    from patankar.layer import layer
    from patankar.food import foodphysics
    from patankar.migration import senspatankar

    fo_grid = make_fo_grid(Fo_max=Fo_max, n_fo=n_fo, fo_min_floor=fo_min_floor)

    # Normalize volume to preserve V_F/V_P ratio
    V_P_real = lP_m * surface_area
    if V_P_real > 0:
        food_volume_normalized = food_volume / V_P_real
    else:
        food_volume_normalized = food_volume

    # Normalize h to preserve Biot number: Bi = h * l / D
    if D_real > 0:
        h_normalized = h * lP_m / D_real
    else:
        h_normalized = h

    # Create normalized layer (l=1, D=1, CP0=1)
    lay = layer(l=1.0, D=1.0, k=k, C0=1.0, T=contact_temperature_degC)

    # Create food physics
    food = foodphysics(
        k=k,
        k0=k0,
        h=h_normalized,
        surfacearea=1.0,
        volume=food_volume_normalized,
        contacttime=float(fo_grid[-1]),
        contacttemperature=contact_temperature_degC,
        CF0=CF0,
    )

    # Solve
    sol = senspatankar(lay, food, t=fo_grid, autotime=False)
    # Route CF extraction through the shared helper (workplan § 1.6).
    # For a full-grid call like this one the solver typically does not
    # extend sol.t, but going through cf_at_user_grid is defence against
    # future call-pattern changes and keeps the survey package uniform.
    from survey.utils import cf_at_user_grid
    cf = cf_at_user_grid(sol, fo_grid)

    return fo_grid, cf


def worker_compute_curve(payload: Dict[str, Any]) -> Dict[str, Any]:
    """
    Worker function for parallel master curve computation.

    Handles cache checking, lock acquisition, and computation.

    Parameters
    ----------
    payload : dict
        Must contain:
        - key_dict: MasterCurveKey as dict
        - cache_dir: Cache directory path

    Returns
    -------
    dict
        Result with status ('hit' or 'miss'), key_hash, fo, cf.
    """
    from survey.cache import MasterCurveCache, MasterCurveKey

    key_dict = payload["key_dict"]
    cache_dir = payload["cache_dir"]
    fo_min_floor = payload.get("fo_min_floor", 1e-15)

    key = MasterCurveKey(**key_dict)
    cache = MasterCurveCache(cache_dir)

    # Check cache
    if cache.exists(key):
        fo, cf = cache.load(key)
        return {
            "status": "hit",
            "key_hash": key.stable_hash(),
            "fo": fo,
            "cf": cf,
        }

    # Try to acquire lock
    acquired = cache.acquire_lock(key)
    if not acquired:
        # Another process is computing; wait for result
        if cache.wait_for_result(key, timeout_s=120.0):
            fo, cf = cache.load(key)
            return {
                "status": "hit",
                "key_hash": key.stable_hash(),
                "fo": fo,
                "cf": cf,
            }
        # Timeout - proceed anyway (better than failing)

    try:
        # Compute master curve
        fo, cf = solve_master_curve(
            k=key.k,
            k0=key.k0,
            lP_m=key.lP_m,
            D_real=key.D,
            Fo_max=key.Fo_max,
            n_fo=key.n_fo,
            contact_temperature_degC=key.contact_temperature_degC,
            h=key.h,
            surface_area=key.surface_area,
            food_volume=key.food_volume,
            CF0=key.CF0,
            fo_min_floor=fo_min_floor,
        )

        # Save to cache
        cache.save(key, fo, cf)

        return {
            "status": "miss",
            "key_hash": key.stable_hash(),
            "fo": fo,
            "cf": cf,
        }
    finally:
        cache.release_lock(key)
