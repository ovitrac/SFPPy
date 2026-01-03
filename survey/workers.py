"""
survey/workers.py — Parallel Worker Functions for Master Curve Computation
==========================================================================

Worker functions compatible with multiprocessing for computing
master curves in parallel.

@project: SFPPy/INSERM — Survey-scale exposure estimation
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

    lo = math.log10(max(fo_min_floor, 1e-30))
    hi = math.log10(max(Fo_max, fo_min_floor))
    grid = np.logspace(lo, hi, int(max(4, n_fo)), dtype=float)
    grid = np.unique(np.concatenate([np.array([0.0], dtype=float), grid]))
    return grid


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
    cf = np.array(sol.CF, dtype=float).reshape(-1)

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
