"""Unified real-units N-layer × K-step migration solver core (WP-1a, 2026-06-22).

One driver for the survey simulation paths. Generalizes the proven real-units bilayer two-step chain
from 2 layers to **N layers**, and from
the hard-wired storage→oven to **1 or 2 steps**. PHYSICAL diffusivities are passed to the kernel; the
engine's own Fourier normalisation does the scaling (C-REALUNITS — the kernel classifies long/short
regime by the magnitude of D it receives, so it must never see D/D_ref). Dimensionless `(fo…, surface)`
are returned ONLY as cache index/metadata; the SOLVE is real-units.

Scope (WP-1a): N spatial layers; K ∈ {1, 2} steps (the only cases in production). Reference layer =
min-permeability `D/(k·l)` at step 1, reused at later steps (same i_ref) — matches the current solvers
exactly so the refactor is numerically identical (A1-NUM). K>2 is intentionally not generalised.

C-RESUME-MESH: step k+1 resumes step k on the SAME mesh (geometry unchanged across steps); topology
changes (layer add/remove) are NOT handled here.

@author: Olivier Vitrac, PhD, HDR — @license: MIT
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Optional, Tuple
import numpy as np


@dataclass(frozen=True)
class LayerSpec:
    """One spatial layer, physical units. D_T2 = None for a single-step run."""
    l_m: float
    D_T1: float
    k: float
    C0: float
    D_T2: Optional[float] = None


def reference_layer(layers: List[LayerSpec]) -> Tuple[int, float, float]:
    """i_ref = highest resistance = min permeability D/(k·l) at step 1 (matches
    _normalise_bilayer / _build_bilayer_tasks). Returns (i_ref, l_ref, D_ref_T1).
    Ties keep the earlier (food-contact) layer, as the current solvers do."""
    best_i, best_perm = 0, float("inf")
    for i, L in enumerate(layers):
        perm = L.D_T1 / (L.k * L.l_m) if L.k * L.l_m > 0 else float("inf")
        if perm < best_perm:
            best_perm, best_i = perm, i
    D_ref_T1 = layers[best_i].D_T1
    if D_ref_T1 <= 0:
        D_ref_T1 = max([L.D_T1 for L in layers] + [1e-30])
    return best_i, layers[best_i].l_m, D_ref_T1


def _build_stack(layers: List[LayerSpec], step: int, T_degC: float):
    """Assemble the N-layer kernel stack with PHYSICAL D for the given step (1 or 2)."""
    from patankar.layer import layer
    stack = None
    for L in layers:
        D = L.D_T1 if step == 1 else L.D_T2
        lay = layer(l=L.l_m, D=D, k=L.k, C0=L.C0, T=T_degC)
        stack = lay if stack is None else stack + lay
    return stack


def solve_chain_master_surface(
    *,
    layers: List[LayerSpec],
    k0: float, h: float, surface_area: float, food_volume: float,
    T1_degC: float, CF0: float,
    Fo1_max: float, n_fo1: int,
    focal_layer: int = 0,
    T2_degC: Optional[float] = None, Fo2_max: Optional[float] = None, n_fo2: Optional[int] = None,
    fo_min_floor: float = 1e-15, eq_rel_tol: float = 1e-3,
):
    """Real-units N-layer master curve (1-step) or surface (2-step).

    Returns:
      1-step  -> (fo1_grid, curve[n1])            curve = CF/CP0 at fo1_grid
      2-step  -> (fo1_grid, fo2_grid, surface[n1,n2])

    Numerically identical to ``solve_twostep_master_surface_v3`` for N=2 (it IS that algorithm with the
    layer count generalised). For N=1 it is the native real-units monolayer chain.
    """
    from patankar.food import foodphysics
    from patankar.migration import senspatankar
    from survey.workers import make_fo_grid
    from survey.workers import _cfeq_from_sol, _eq_knee
    from survey.utils import cf_at_user_grid

    two_step = T2_degC is not None
    if two_step and (Fo2_max is None or n_fo2 is None):
        raise ValueError("two-step requires Fo2_max and n_fo2")
    for L in layers:
        if not np.isfinite(L.D_T1) or L.D_T1 <= 0 or (two_step and (L.D_T2 is None or L.D_T2 <= 0)):
            raise ValueError("C-UNITS: every layer needs real physical D (>0); never D/D_ref")

    i_ref, l_ref, D_ref_T1 = reference_layer(layers)
    food_k = layers[focal_layer].k

    fo1_grid = make_fo_grid(Fo_max=Fo1_max, n_fo=n_fo1, fo_min_floor=fo_min_floor)
    t1_phys = fo1_grid * (l_ref ** 2 / D_ref_T1)        # dimensionless axis → physical time
    stack_T1 = _build_stack(layers, 1, T1_degC)

    # ---------- single-step ----------
    if not two_step:
        t1_solve = t1_phys[t1_phys > 0.0]
        food = foodphysics(k=food_k, k0=k0, h=(h, "m/s"), surfacearea=(surface_area, "m**2"),
                           volume=(food_volume, "m**3"), contacttime=(float(t1_solve[-1]), "s"),
                           contacttemperature=(T1_degC, "degC"), CF0=CF0)
        sol = senspatankar(stack_T1, food, t=t1_solve, autotime=False, timescale="sqrt")
        cf = np.asarray(cf_at_user_grid(sol, t1_solve), dtype=float)
        cfeq = _cfeq_from_sol(sol)
        knee = _eq_knee(cf, cfeq, eq_rel_tol)
        if cfeq > 0 and knee < len(cf):
            cf[knee:] = cfeq
        curve = np.empty(len(fo1_grid))
        curve[0] = cf[0] if fo1_grid[0] > 0 else 0.0
        curve[1:] = cf[:len(fo1_grid) - 1]
        return fo1_grid, curve

    # ---------- two-step ----------
    D_ref_T2 = layers[i_ref].D_T2
    if D_ref_T2 <= 0:
        D_ref_T2 = max([L.D_T2 for L in layers] + [1e-30])
    fo2_grid = make_fo_grid(Fo_max=Fo2_max, n_fo=n_fo2, fo_min_floor=fo_min_floor)
    t2_phys = fo2_grid * (l_ref ** 2 / D_ref_T2)
    t2_solve = t2_phys[t2_phys > 0.0]              # fo2=0 is the storage-alone carry column

    stack_T2 = _build_stack(layers, 2, T2_degC)
    food_T2 = foodphysics(k=food_k, k0=k0, h=(h, "m/s"), surfacearea=(surface_area, "m**2"),
                          volume=(food_volume, "m**3"), contacttime=(float(t2_solve[-1]), "s"),
                          contacttemperature=(T2_degC, "degC"), CF0=CF0)

    surface = np.empty((len(fo1_grid), len(fo2_grid)), dtype=float)
    cfeq_T1 = cfeq_T2 = None
    sol1_frozen = cf1_frozen = None

    for i, t1 in enumerate(t1_phys):
        t1_use = float(t1) if t1 > 0 else float(t2_solve[0])
        if sol1_frozen is not None:
            sol1, cf1 = sol1_frozen, cf1_frozen
        else:
            food_T1 = foodphysics(k=food_k, k0=k0, h=(h, "m/s"), surfacearea=(surface_area, "m**2"),
                                  volume=(food_volume, "m**3"), contacttime=(t1_use, "s"),
                                  contacttemperature=(T1_degC, "degC"), CF0=CF0)
            sol1 = senspatankar(stack_T1, food_T1, t=np.array([0.0, t1_use]),
                                autotime=False, timescale="sqrt")
            cf1 = float(sol1.CFtarget)
            if cfeq_T1 is None:
                cfeq_T1 = _cfeq_from_sol(sol1)
            if cfeq_T1 > 0 and abs(cf1 - cfeq_T1) <= eq_rel_tol * cfeq_T1:
                sol1_frozen, cf1_frozen = sol1, cf1
        surface[i, 0] = cf1                                    # fo2=0 → storage alone
        sol2 = sol1.resume(multilayer=stack_T2, medium=food_T2, CF0=cf1,
                           t=t2_solve, autotime=False)
        cf_row = np.asarray(cf_at_user_grid(sol2, t2_solve), dtype=float)
        if cfeq_T2 is None:
            cfeq_T2 = _cfeq_from_sol(sol2)
        if cfeq_T2 and cfeq_T2 > 0:
            knee = _eq_knee(cf_row, cfeq_T2, eq_rel_tol)
            if knee < len(cf_row):
                cf_row[knee:] = cfeq_T2
        surface[i, 1:] = cf_row[:surface.shape[1] - 1]

    return fo1_grid, fo2_grid, surface


_T2_FLOOR_S = 1e-3   # tiny physical time standing in for "oven ~ off" (mirrors mono_twostep_chained)


def chain_at_times(
    *,
    layers: List[LayerSpec],
    k0: float, h: float, surface_area: float, food_volume: float,
    T1_degC: float, CF0: float,
    t1_grid,
    focal_layer: int = 0,
    T2_degC: Optional[float] = None, t2_grid=None,
    layers_step2: Optional[List[LayerSpec]] = None,
    eq_rel_tol: float = 1e-3,
    _legacy_twostep: bool = False,
    c0eq_out: Optional[dict] = None,
):
    """PHYSICAL chain evaluated at the scenario's ACTUAL time-prior points — no Fo grid, no
    interpolation (the physical-primary query primitive, WP-physical P1). N-layer generalisation of
    ``mono_twostep_chained.chain_substance``; PHYSICAL D throughout.

    Returns (C0=1 → CF == g = CF/CP0):
      1-step  -> g1[n1]           = CF at each t1 in t1_grid
      2-step  -> (g1[n1], g12[n1,n2])

    2-step mirrors chain_substance exactly: storage solve per t1 (timescale='sqrt'), then one
    resume over the oven grid (multi-t2), re-passing CF0=sol1.CFtarget (food carry), with the oven
    DURATION integrated from ~0 (``_T2_FLOOR_S``) so the 0→t2 transient is captured and CF is read at
    the real t2 durations. CF12 ≥ CF1 by construction.

    OVERPACK REMOVAL (WP-6, OV 2026-06-25): ``layers_step2`` overrides the layer set used for the OVEN
    step (step 1 always uses the full ``layers`` — the overpack is present during storage). It is the
    RETAINED layer set after dropping any overpack paper/board layer(s):
      - ``None``  → no removal; the oven uses the full ``layers`` (current behaviour).
      - non-empty → the oven solves on the reduced stack; the storage profile is remapped onto the
        shorter mesh by ``senspatankar`` (``Cxprevious.interp`` over the new, shorter geometry), so the
        removed OUTER layer's absorbed mass is discarded — it leaves with the carton, not redistributed.
        Requires the removed layer to be OUTER (high-x); the food-contact/focal layer must be retained.
      - empty list → monolayer overpack: nothing is left in the oven, so there is no further migration —
        the food stays at its storage value (g12[:, :] = g1[:, None]).
    """
    from patankar.food import foodphysics
    from patankar.migration import senspatankar
    from survey.workers import _cfeq_from_sol, _eq_knee
    from survey.utils import cf_at_user_grid

    two_step = T2_degC is not None
    for L in layers:
        if not np.isfinite(L.D_T1) or L.D_T1 <= 0 or (two_step and (L.D_T2 is None or L.D_T2 <= 0)):
            raise ValueError("C-UNITS: every layer needs real physical D (>0); never D/D_ref")
    food_k = layers[focal_layer].k
    t1_arr = np.asarray(t1_grid, dtype=float)
    stack_T1 = _build_stack(layers, 1, T1_degC)

    # ---------- single-step: one multi-t solve over the real t1 grid ----------
    if not two_step:
        t1_solve = t1_arr[t1_arr > 0.0]
        food = foodphysics(k=food_k, k0=k0, h=(h, "m/s"), surfacearea=(surface_area, "m**2"),
                           volume=(food_volume, "m**3"), contacttime=(float(t1_solve[-1]), "s"),
                           contacttemperature=(T1_degC, "degC"), CF0=CF0)
        sol = senspatankar(stack_T1, food, t=t1_solve, autotime=False, timescale="sqrt")
        if c0eq_out is not None:                       # 1-step: genuine equilibrium scale
            c0eq_out["C0eq"] = float(np.ravel(sol.C0eq)[0])
            c0eq_out["C0eq_is_default"] = bool(getattr(sol, "C0eq_is_default", False))
        cf = np.asarray(cf_at_user_grid(sol, t1_arr), dtype=float)
        cfeq = _cfeq_from_sol(sol)
        if cfeq > 0:
            knee = _eq_knee(cf, cfeq, eq_rel_tol)
            if knee < len(cf):
                cf[knee:] = cfeq
        return cf

    # ---------- two-step: per-t1 storage + multi-t2 resume (== chain_substance) ----------
    t2_arr = np.asarray(t2_grid, dtype=float)
    _pts = t2_arr[t2_arr > _T2_FLOOR_S]
    t2_solve = np.concatenate(([_T2_FLOOR_S], _pts))     # strictly increasing, starts at ~0
    n1, n2 = len(t1_arr), len(t2_arr)
    g1 = np.empty(n1); g12 = np.empty((n1, n2))

    # Overpack removal: the OVEN step uses the retained layer set (step 1 keeps the full stack).
    oven_layers = layers if layers_step2 is None else layers_step2
    overpack_only = (layers_step2 is not None and len(layers_step2) == 0)
    if not overpack_only:
        stack_T2 = _build_stack(oven_layers, 2, T2_degC)
        food_T2 = foodphysics(k=food_k, k0=k0, h=(h, "m/s"), surfacearea=(surface_area, "m**2"),
                              volume=(food_volume, "m**3"), contacttime=(float(t2_solve[-1]), "s"),
                              contacttemperature=(T2_degC, "degC"), CF0=CF0)

    if _legacy_twostep:
        # ---- LEGACY path (kept for falsification): one storage solve PER t1 ----
        for i, t1 in enumerate(t1_arr):
            t1_use = float(t1) if t1 > 0 else _T2_FLOOR_S
            food_T1 = foodphysics(k=food_k, k0=k0, h=(h, "m/s"), surfacearea=(surface_area, "m**2"),
                                  volume=(food_volume, "m**3"), contacttime=(t1_use, "s"),
                                  contacttemperature=(T1_degC, "degC"), CF0=CF0)
            s1 = senspatankar(stack_T1, food_T1, t=np.array([0.0, t1_use]),
                              autotime=False, timescale="sqrt")
            g1[i] = float(s1.CFtarget)
            if overpack_only:
                g12[i, :] = g1[i]
                continue
            s2 = s1.resume(multilayer=stack_T2, medium=food_T2, CF0=float(s1.CFtarget),
                           t=t2_solve, autotime=False)
            g12[i, :] = np.asarray(cf_at_user_grid(s2, t2_arr), dtype=float)[:n2]
        return g1, g12

    # ---- OPTIMISED path: per-t1 storage solve (accurate — identical to the legacy
    # storage step) + TAIL-INTEGRATED oven resume. The 17x win comes entirely from
    # the oven tail-integration (exact, 0 error vs no-tail); the oven was the whole
    # cost. The earlier storage single-march was reverted: it added a small error on
    # sharp frozen ICs for ~0 gain (storage is already cheap). resumeat(t1_use, ...)
    # restarts from the storage end-state (= resume) but activates tail_integration.
    for i, t1 in enumerate(t1_arr):
        t1_use = float(t1) if t1 > 0 else _T2_FLOOR_S
        food_T1 = foodphysics(k=food_k, k0=k0, h=(h, "m/s"), surfacearea=(surface_area, "m**2"),
                              volume=(food_volume, "m**3"), contacttime=(t1_use, "s"),
                              contacttemperature=(T1_degC, "degC"), CF0=CF0)
        s1 = senspatankar(stack_T1, food_T1, t=np.array([0.0, t1_use]),
                          autotime=False, timescale="sqrt")
        if c0eq_out is not None and "C0eq" not in c0eq_out:   # capture once (same across t1)
            # record the genuine STORAGE equilibrium scale; flag True because the
            # chain resumes the oven (a defaulted-C0eq step) — C0eq is therefore
            # not a single closed-system equilibrium for the whole two-step chain.
            c0eq_out["C0eq"] = float(np.ravel(s1.C0eq)[0])
            c0eq_out["C0eq_is_default"] = True
        g1[i] = float(s1.CFtarget)
        if overpack_only:
            g12[i, :] = g1[i]
            continue
        # tail_integration=True: opt-in fast hot-step (relaxed late-Fo tail).
        # senspatankar's own default stays legacy (False); the survey activates it.
        s2 = s1.resumeat(t1_use, multilayer=stack_T2, medium=food_T2,
                         t=t2_solve, autotime=False, tail_integration=True)
        g12[i, :] = np.asarray(cf_at_user_grid(s2, t2_arr), dtype=float)[:n2]
    return g1, g12
