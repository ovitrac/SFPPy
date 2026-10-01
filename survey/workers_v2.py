"""
survey/workers_v2.py — Bilayer and Two-Step Master Curve Computation
=====================================================================

Extends workers.py for:
- Bilayer packaging (2 layers with different D, k, l, C0)
- Two-step contact (chained simulations at different temperatures)
- Combined bilayer + two-step

The monolayer single-step path is untouched (still in workers.py).
Bilayer and two-step results use a separate cache subdirectory.

@project: SFPPy — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

import os
import json
import hashlib
import math
from pathlib import Path
from typing import Dict, Any, Tuple, List, Optional
from dataclasses import dataclass, asdict

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

from survey.workers import make_fo_grid, assert_bilayer_mass_conservation


# =====================================================================
# Bilayer Master Curve
# =====================================================================

@dataclass(frozen=True)
class BilayerCurveKey:
    """
    Key that uniquely identifies a bilayer master curve.

    Includes both layers' properties and the focal layer index.
    """
    # Layer 1 (food contact)
    polymer_1: str
    D_1: float
    k_1: float
    l_1_m: float
    C0_1: float      # 1.0 for focal layer, 0.0 for non-focal

    # Layer 2 (outer/barrier)
    polymer_2: str
    D_2: float
    k_2: float
    l_2_m: float
    C0_2: float

    # Shared
    k0: float
    h: float
    surface_area: float
    food_volume: float
    contact_temperature_degC: float
    CF0: float
    Fo_max: float
    n_fo: int
    focal_layer: int   # 0 or 1: which layer's CP0 is the variable
    # 'v2.1': food-volume normalisation fix (V_ref = l_ref·A). Invalidates
    # mass-non-conserving curves cached under the unversioned key.
    worker_version: str = 'v2.1'

    def stable_hash(self) -> str:
        """Deterministic hash for cache."""
        d = asdict(self)
        raw = json.dumps(d, sort_keys=True, default=str)
        return hashlib.sha256(raw.encode()).hexdigest()[:16]


def solve_bilayer_master_curve(
    *,
    # Layer 1 (food contact side)
    polymer_1: str, D_1: float, k_1: float, l_1_m: float, C0_1: float,
    # Layer 2 (outer/barrier)
    polymer_2: str, D_2: float, k_2: float, l_2_m: float, C0_2: float,
    # Shared physics
    k0: float, h: float, surface_area: float, food_volume: float,
    contact_temperature_degC: float, CF0: float,
    Fo_max: float, n_fo: int,
    focal_layer: int = 0,
    fo_min_floor: float = 1e-15,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Solve bilayer master curve g(Fo) = CF / CP0_focal.

    Uses the reference layer (highest resistance = min D/(k*l))
    to define the Fo timescale. Both layers are built with their
    real (normalized) properties.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (fo_grid, cf_over_cp0) where cf_over_cp0 = CF / CP0_focal
    """
    from patankar.layer import layer
    from patankar.food import foodphysics
    from patankar.migration import senspatankar

    # Reference layer: the one with lowest D/(k*l) = highest resistance
    perm_1 = D_1 / (k_1 * l_1_m) if k_1 * l_1_m > 0 else float('inf')
    perm_2 = D_2 / (k_2 * l_2_m) if k_2 * l_2_m > 0 else float('inf')
    if perm_1 <= perm_2:
        i_ref = 0
        l_ref, D_ref = l_1_m, D_1
    else:
        i_ref = 1
        l_ref, D_ref = l_2_m, D_2

    if D_ref <= 0:
        D_ref = max(D_1, D_2, 1e-30)

    # Normalize relative to reference layer
    l_1_norm = l_1_m / l_ref
    l_2_norm = l_2_m / l_ref
    D_1_norm = D_1 / D_ref
    D_2_norm = D_2 / D_ref

    # Volume ratio preservation — normalise food on the SAME l_ref·A scale as
    # the layers (NOT total packaging volume), so the kernel dilution number
    # L = A·l_ref/V_F (migration.py:4168) is preserved. Using (l_1+l_2)·A
    # inflated CF/CP0 by (l_1+l_2)/l_ref and broke mass conservation.
    V_ref = l_ref * surface_area
    food_volume_norm = food_volume / V_ref if V_ref > 0 else food_volume

    # Biot number preservation (using reference layer)
    h_norm = h * l_ref / D_ref if D_ref > 0 else h

    # Fo grid (reference layer timescale)
    fo_grid = make_fo_grid(Fo_max=Fo_max, n_fo=n_fo, fo_min_floor=fo_min_floor)

    T = contact_temperature_degC

    # Build bilayer: layer 1 (food contact) + layer 2 (outer)
    lay1 = layer(l=l_1_norm, D=D_1_norm, k=k_1, C0=C0_1, T=T)
    lay2 = layer(l=l_2_norm, D=D_2_norm, k=k_2, C0=C0_2, T=T)
    bilayer = lay1 + lay2

    # Food
    food = foodphysics(
        k=k_1 if focal_layer == 0 else k_2,
        k0=k0,
        h=h_norm,
        surfacearea=1.0,
        volume=food_volume_norm,
        contacttime=float(fo_grid[-1]),
        contacttemperature=T,
        CF0=CF0,
    )

    sol = senspatankar(bilayer, food, t=fo_grid, autotime=False)
    cf = np.array(sol.CF, dtype=float).reshape(-1)

    assert_bilayer_mass_conservation(
        cf, l_1_m=l_1_m, l_2_m=l_2_m, C0_1=C0_1, C0_2=C0_2,
        surface_area=surface_area, food_volume=food_volume, context="(1-step v2)",
    )
    return fo_grid, cf


# =====================================================================
# Two-Step Master Surface (2D)
# =====================================================================

@dataclass(frozen=True)
class TwoStepSurfaceKey:
    """
    Key for a 2D master surface g(Fo1, Fo2) for two-step contact.
    """
    # Layer stack (bilayer or monolayer)
    polymer_1: str
    D_1_T1: float       # D at T1
    D_1_T2: float       # D at T2
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
    # 'v2.1': food-volume normalisation fix (V_ref = l_ref·A). Invalidates
    # mass-non-conserving surfaces cached under the unversioned key.
    worker_version: str = 'v2.1'

    def stable_hash(self) -> str:
        d = asdict(self)
        raw = json.dumps(d, sort_keys=True, default=str)
        return hashlib.sha256(raw.encode()).hexdigest()[:16]


def solve_twostep_master_surface(
    *,
    # Layer 1
    polymer_1: str, D_1_T1: float, D_1_T2: float,
    k_1: float, l_1_m: float, C0_1: float,
    # Layer 2
    polymer_2: str, D_2_T1: float, D_2_T2: float,
    k_2: float, l_2_m: float, C0_2: float,
    # Shared
    k0: float, h: float, surface_area: float, food_volume: float,
    T1_degC: float, T2_degC: float, CF0: float,
    Fo1_max: float, Fo2_max: float,
    n_fo1: int = 30, n_fo2: int = 15,
    focal_layer: int = 0,
    fo_min_floor: float = 1e-15,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute 2D master surface g(Fo1, Fo2) for two-step contact.

    Step 1 at T1 for Fo1 range, then step 2 at T2 for Fo2 range.
    Uses >> chaining to propagate CF0 from step 1 to step 2.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        (fo1_grid, fo2_grid, surface) where surface[i,j] = CF/CP0
        at (Fo1=fo1_grid[i], Fo2=fo2_grid[j])
    """
    from patankar.layer import layer
    from patankar.food import foodphysics
    from patankar.migration import senspatankar

    # Reference layer at T1 for Fo1 timescale
    perm_1 = D_1_T1 / (k_1 * l_1_m) if k_1 * l_1_m > 0 else float('inf')
    perm_2 = D_2_T1 / (k_2 * l_2_m) if k_2 * l_2_m > 0 else float('inf')
    if perm_1 <= perm_2:
        l_ref, D_ref_T1, D_ref_T2 = l_1_m, D_1_T1, D_1_T2
    else:
        l_ref, D_ref_T1, D_ref_T2 = l_2_m, D_2_T1, D_2_T2

    if D_ref_T1 <= 0:
        D_ref_T1 = max(D_1_T1, D_2_T1, 1e-30)
    if D_ref_T2 <= 0:
        D_ref_T2 = max(D_1_T2, D_2_T2, 1e-30)

    # Normalize at T1
    l_1_norm = l_1_m / l_ref
    l_2_norm = l_2_m / l_ref

    # Normalise food on the l_ref·A scale (same as the layers), preserving the
    # kernel dilution number L = A·l_ref/V_F. (Was (l_1+l_2)·A — see fix above.)
    V_ref = l_ref * surface_area
    VF_norm = food_volume / V_ref if V_ref > 0 else food_volume
    h_norm_T1 = h * l_ref / D_ref_T1 if D_ref_T1 > 0 else h

    fo1_grid = make_fo_grid(Fo1_max, n_fo1, fo_min_floor)
    fo2_grid = make_fo_grid(Fo2_max, n_fo2, fo_min_floor)

    surface = np.zeros((len(fo1_grid), len(fo2_grid)), dtype=float)

    for i, Fo1 in enumerate(fo1_grid):
        # Build bilayer at T1
        D_1n_T1 = D_1_T1 / D_ref_T1
        D_2n_T1 = D_2_T1 / D_ref_T1
        lay1_T1 = layer(l=l_1_norm, D=D_1n_T1, k=k_1, C0=C0_1, T=T1_degC)
        lay2_T1 = layer(l=l_2_norm, D=D_2n_T1, k=k_2, C0=C0_2, T=T1_degC)
        bi_T1 = lay1_T1 + lay2_T1

        food1 = foodphysics(
            k=k_1 if focal_layer == 0 else k_2,
            k0=k0, h=h_norm_T1, surfacearea=1.0, volume=VF_norm,
            contacttime=float(max(Fo1, 1e-20)),
            contacttemperature=T1_degC, CF0=CF0,
        )

        if Fo1 <= 0:
            # No step 1 — go directly to step 2
            CF0_step2 = CF0
            Cx_prev = None
        else:
            sol1 = senspatankar(bi_T1, food1, t=np.array([0.0, float(Fo1)]),
                                autotime=False)
            # sol.CFtarget is identical to sol.CF[-1] here because
            # food1.contacttime == max(t) == Fo1 (no post-contact window
            # extension in this call pattern). CFtarget is the safer
            # accessor by convention — see survey/utils/cf_extract.py.
            CF0_step2 = float(sol1.CFtarget)
            Cx_prev = None  # simplified: use CF0 propagation only

        # Step 2: bilayer at T2
        h_norm_T2 = h * l_ref / D_ref_T2 if D_ref_T2 > 0 else h
        D_1n_T2 = D_1_T2 / D_ref_T2
        D_2n_T2 = D_2_T2 / D_ref_T2

        # After step 1, polymer is partially depleted. For the normalized
        # master surface, we approximate: C0 in polymer remains ~C0_initial
        # (conservative, since depletion is small for most survey scenarios).
        lay1_T2 = layer(l=l_1_norm, D=D_1n_T2, k=k_1, C0=C0_1, T=T2_degC)
        lay2_T2 = layer(l=l_2_norm, D=D_2n_T2, k=k_2, C0=C0_2, T=T2_degC)
        bi_T2 = lay1_T2 + lay2_T2

        for j, Fo2 in enumerate(fo2_grid):
            if Fo2 <= 0:
                surface[i, j] = CF0_step2
                continue

            # Rescale Fo2 to T2 reference timescale
            # Fo2 is in T1 reference units; convert to T2 solver time
            # t_real = Fo2 * l_ref² / D_ref_T1 (Fo2 defined via T1 ref)
            # solver_time_T2 = t_real / (l_ref² / D_ref_T2) = Fo2 * D_ref_T2/D_ref_T1
            Fo2_T2 = Fo2 * D_ref_T2 / D_ref_T1 if D_ref_T1 > 0 else Fo2

            food2 = foodphysics(
                k=k_1 if focal_layer == 0 else k_2,
                k0=k0, h=h_norm_T2, surfacearea=1.0, volume=VF_norm,
                contacttime=float(max(Fo2_T2, 1e-20)),
                contacttemperature=T2_degC, CF0=CF0_step2,
            )

            sol2 = senspatankar(bi_T2, food2,
                                t=np.array([0.0, float(Fo2_T2)]),
                                autotime=False)
            # Same convention as CF0_step2 above — contacttime == max(t),
            # so CFtarget == CF[-1] here; we prefer CFtarget.
            surface[i, j] = float(sol2.CFtarget)

    assert_bilayer_mass_conservation(
        surface, l_1_m=l_1_m, l_2_m=l_2_m, C0_1=C0_1, C0_2=C0_2,
        surface_area=surface_area, food_volume=food_volume, context="(2-step v2)",
    )
    return fo1_grid, fo2_grid, surface


# =====================================================================
# Cache for Bilayer and Two-Step
# =====================================================================

class BilayerCache:
    """
    Persistent cache for bilayer master curves.

    Separate from the monolayer MasterCurveCache to avoid collisions.
    Layout: cache_dir/bilayer_curves/<hash>.npz
    """

    def __init__(self, cache_dir: str):
        self.root = Path(cache_dir).expanduser().resolve() / "bilayer_curves"
        self.root.mkdir(parents=True, exist_ok=True)
        self._stats = {"hits": 0, "misses": 0}

    def _path(self, key) -> Path:
        return self.root / f"{key.stable_hash()}.npz"

    def exists(self, key) -> bool:
        return self._path(key).exists()

    def load_1d(self, key: BilayerCurveKey) -> Tuple[np.ndarray, np.ndarray]:
        data = np.load(self._path(key))
        self._stats["hits"] += 1
        return data["fo_grid"], data["cf_over_cp0"]

    def save_1d(self, key: BilayerCurveKey, fo: np.ndarray, cf: np.ndarray):
        np.savez_compressed(self._path(key), fo_grid=fo, cf_over_cp0=cf)
        meta = self.root / f"{key.stable_hash()}.json"
        meta.write_text(json.dumps(asdict(key), indent=2, sort_keys=True, default=str))
        self._stats["misses"] += 1

    def load_2d(self, key: TwoStepSurfaceKey) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        data = np.load(self._path(key))
        self._stats["hits"] += 1
        return data["fo1_grid"], data["fo2_grid"], data["surface"]

    def save_2d(self, key: TwoStepSurfaceKey, fo1: np.ndarray, fo2: np.ndarray,
                surface: np.ndarray):
        np.savez_compressed(self._path(key), fo1_grid=fo1, fo2_grid=fo2, surface=surface)
        meta = self.root / f"{key.stable_hash()}.json"
        meta.write_text(json.dumps(asdict(key), indent=2, sort_keys=True, default=str))
        self._stats["misses"] += 1

    @property
    def stats(self) -> Dict[str, int]:
        return dict(self._stats)


# =====================================================================
# Worker Entry Points (for parallel processing)
# =====================================================================

def worker_bilayer_curve(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Worker for parallel bilayer master curve computation with caching."""
    key = BilayerCurveKey(**payload["key_dict"])
    cache = BilayerCache(payload["cache_dir"])

    if cache.exists(key):
        fo, cf = cache.load_1d(key)
        return {"status": "hit", "key_hash": key.stable_hash(), "fo": fo, "cf": cf}

    fo, cf = solve_bilayer_master_curve(**{
        k: v for k, v in asdict(key).items()
        if k not in ('focal_layer', 'worker_version')
    }, focal_layer=key.focal_layer)

    cache.save_1d(key, fo, cf)
    return {"status": "miss", "key_hash": key.stable_hash(), "fo": fo, "cf": cf}


def worker_twostep_surface(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Worker for parallel two-step surface computation with caching."""
    key = TwoStepSurfaceKey(**payload["key_dict"])
    cache = BilayerCache(payload["cache_dir"])

    if cache.exists(key):
        fo1, fo2, surf = cache.load_2d(key)
        return {"status": "hit", "key_hash": key.stable_hash(),
                "fo1": fo1, "fo2": fo2, "surface": surf}

    fo1, fo2, surf = solve_twostep_master_surface(**{
        k: v for k, v in asdict(key).items()
        if k not in ('focal_layer', 'worker_version')
    }, focal_layer=key.focal_layer)

    cache.save_2d(key, fo1, fo2, surf)
    return {"status": "miss", "key_hash": key.stable_hash(),
            "fo1": fo1, "fo2": fo2, "surface": surf}
