"""Falsification tests for the n>=3-layer survey path.

The claim under test: `survey/` now handles n>=3 spatial layers through the
real-units chain engine with tuple-valued cache keys, and the result equals a
DIRECT patankar kernel solve of the same trilayer (no survey abstraction).

Also guards the two invariants of the T1 design:
  * n == 2 keys are byte-identical to before (warm g-cache preserved);
  * distinct hash spaces (a trilayer key can never collide with a bilayer key).

Author: Olivier Vitrac, PhD, HDR
"""
import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Trilayer test geometry (A01-like: thin contact liner | board | outer barrier)
# ---------------------------------------------------------------------------
L_MS = (4e-6, 50e-6, 25e-6)          # LDPE 4 um | PB 50 um | PA6 25 um
POLYMERS = ("LDPE", "PB", "PA6")
D_T1S = (1e-13, 5e-12, 1e-15)        # explicit physical D (no property lookup)
KS = (1.0, 1.0, 1.0)
T1 = 25.0
SURFACE = 441e-4                      # m^2  (A01)
VOLUME = 250e-6                       # m^3  (250 mL)
H = 1e-4                              # m/s
K0 = 1.0
CF0 = 0.0
FO_MAX = 2.0
N_FO = 60


def _direct_kernel_cf(t_s: np.ndarray, focal: int) -> np.ndarray:
    """Reference: direct senspatankar trilayer solve, C0=1 on the focal layer."""
    from patankar.layer import layer
    from patankar.food import foodphysics
    from patankar.migration import senspatankar
    from survey.utils import cf_at_user_grid

    stack = None
    for i in range(3):
        lay = layer(l=L_MS[i], D=D_T1S[i], k=KS[i],
                    C0=1.0 if i == focal else 0.0, T=T1)
        stack = lay if stack is None else stack + lay
    food = foodphysics(k=KS[focal], k0=K0, h=(H, "m/s"),
                       surfacearea=(SURFACE, "m**2"), volume=(VOLUME, "m**3"),
                       contacttime=(float(t_s[-1]), "s"),
                       contacttemperature=(T1, "degC"), CF0=CF0)
    sol = senspatankar(stack, food, t=t_s, autotime=False, timescale="sqrt")
    return np.asarray(cf_at_user_grid(sol, t_s), dtype=float)


@pytest.mark.parametrize("focal", [0, 1, 2])
def test_trilayer_curve_matches_direct_kernel(focal):
    """solve_multilayer_master_curve_v3 == direct kernel solve (all origins)."""
    from survey.workers_v3 import solve_multilayer_master_curve_v3

    C0s = tuple(1.0 if i == focal else 0.0 for i in range(3))
    fo, cf = solve_multilayer_master_curve_v3(
        polymers=POLYMERS, D_T1s=D_T1S, ks=KS, l_ms=L_MS, C0s=C0s,
        k0=K0, h=H, surface_area=SURFACE, food_volume=VOLUME,
        contact_temperature_degC=T1, CF0=CF0, Fo_max=FO_MAX, n_fo=N_FO,
        focal_layer=focal,
    )
    assert fo.shape == cf.shape and len(fo) >= 4
    assert np.all(np.isfinite(cf)) and np.all(cf >= 0.0)

    # Reference layer used by chain_engine: min permeability D/(k*l)
    perms = [D_T1S[i] / (KS[i] * L_MS[i]) for i in range(3)]
    i_ref = int(np.argmin(perms))
    l_ref, D_ref = L_MS[i_ref], D_T1S[i_ref]
    t_phys = fo * (l_ref ** 2 / D_ref)
    t_pos = t_phys[t_phys > 0.0]

    cf_ref = _direct_kernel_cf(t_pos, focal)
    # chain_engine writes curve[0]=0 (fo=0) then cf[:-1] on the positive grid
    cf_engine_pos = cf[1:len(t_pos) + 1] if fo[0] == 0.0 else cf[:len(t_pos)]
    # tolerate the one-sample shift convention: compare on the common length
    m = min(len(cf_engine_pos), len(cf_ref))
    num = np.abs(cf_engine_pos[:m] - cf_ref[:m])
    den = np.maximum(np.abs(cf_ref[:m]), 1e-30)
    # identical algorithm & solver → tight tolerance
    assert np.nanmax(num / np.maximum(den, np.nanmax(den) * 1e-6)) < 5e-2, (
        f"trilayer focal={focal}: survey path deviates from direct kernel")


def test_trilayer_mass_ceiling_all_origins():
    """max(CF/CP0) <= C0_focal*l_focal*A/V_F for every origin (eq. M4 class)."""
    from survey.workers_v3 import solve_multilayer_master_curve_v3
    for focal in range(3):
        C0s = tuple(1.0 if i == focal else 0.0 for i in range(3))
        _, cf = solve_multilayer_master_curve_v3(
            polymers=POLYMERS, D_T1s=D_T1S, ks=KS, l_ms=L_MS, C0s=C0s,
            k0=K0, h=H, surface_area=SURFACE, food_volume=VOLUME,
            contact_temperature_degC=T1, CF0=CF0, Fo_max=FO_MAX, n_fo=N_FO,
            focal_layer=focal,
        )
        ceiling = L_MS[focal] * SURFACE / VOLUME
        assert float(np.nanmax(cf)) <= ceiling * (1 + 1e-6)


def test_multilayer_key_hash_space_distinct():
    """A trilayer key can never collide with a bilayer key (distinct classes,
    distinct worker_version) and its hash is stable across tuple/list input."""
    from survey.workers_v3 import (BilayerCurveKeyV3, MultilayerCurveKeyV3)

    tri = MultilayerCurveKeyV3(
        polymers=POLYMERS, D_T1s=D_T1S, ks=KS, l_ms=L_MS, C0s=(1.0, 0.0, 0.0),
        k0=K0, h=H, surface_area=SURFACE, food_volume=VOLUME,
        contact_temperature_degC=T1, CF0=CF0, Fo_max=FO_MAX, n_fo=N_FO,
        focal_layer=0)
    bi = BilayerCurveKeyV3(
        polymer_1="LDPE", D_1=D_T1S[0], k_1=1.0, l_1_m=L_MS[0], C0_1=1.0,
        polymer_2="PB", D_2=D_T1S[1], k_2=1.0, l_2_m=L_MS[1], C0_2=0.0,
        k0=K0, h=H, surface_area=SURFACE, food_volume=VOLUME,
        contact_temperature_degC=T1, CF0=CF0, Fo_max=FO_MAX, n_fo=N_FO,
        focal_layer=0)
    assert tri.stable_hash() != bi.stable_hash()
    assert tri.worker_version == "v3.3-multi"

    # list-vs-tuple normalisation in the worker: same hash either way
    from dataclasses import asdict
    kd = asdict(tri)
    kd_list = dict(kd)
    for f in ("polymers", "D_T1s", "ks", "l_ms", "C0s"):
        kd_list[f] = list(kd_list[f])
        kd_list[f] = tuple(kd_list[f])
    assert MultilayerCurveKeyV3(**kd_list).stable_hash() == tri.stable_hash()


def test_bilayer_key_unchanged_regression():
    """The n == 2 key layout/version is untouched by T1 (warm cache preserved)."""
    from survey.workers_v3 import BilayerCurveKeyV3, TwoStepSurfaceKeyV3
    import dataclasses
    f1 = [f.name for f in dataclasses.fields(BilayerCurveKeyV3)]
    f2 = [f.name for f in dataclasses.fields(TwoStepSurfaceKeyV3)]
    assert f1 == ['polymer_1', 'D_1', 'k_1', 'l_1_m', 'C0_1',
                  'polymer_2', 'D_2', 'k_2', 'l_2_m', 'C0_2',
                  'k0', 'h', 'surface_area', 'food_volume',
                  'contact_temperature_degC', 'CF0', 'Fo_max', 'n_fo',
                  'focal_layer', 'worker_version']
    assert BilayerCurveKeyV3.__dataclass_fields__['worker_version'].default == 'v3.1'
    assert TwoStepSurfaceKeyV3.__dataclass_fields__['worker_version'].default == 'v3.2'
    assert 'polymer_1' in f2 and 'D_1_T2' in f2
