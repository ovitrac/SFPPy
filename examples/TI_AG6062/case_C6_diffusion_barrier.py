"""
===============================================================================
Case C6 — Barrière fonctionnelle à effet diffusion (tier M3 + FB)
===============================================================================
Case study C6 (section 4.6) of:

    Vitrac O., Nguyen P.-M. Conception sûre des emballages – Cas d'étude et
    outils numériques. Techniques de l'Ingénieur, AG 6 062 (2027, à paraître).

Authors:
    Olivier Vitrac       Innovation Lab, Adservio Group, Puteaux, France
    Phuong-Mai Nguyen    LNE (Laboratoire national de métrologie et d'essais),
                         Trappes, France

Series "Conception sûre des emballages" (Techniques de l'Ingénieur):
    1/3  AG 6 058 (2026) – Principes, modélisation et fondements
    2/3  AG 6 059 – Tiers avancés M2–M3 : partage, diffusion et barrières
    3/3  AG 6 062 – Cas d'étude et outils numériques

Recycled PP (300 µm) contaminated with toluene (C0 = 200 mg/kg), 10 d at 40 °C,
SML = 0.6 mg/kg.
A plasticized PET (wPET) functional barrier on the food-contact side delays
migration. Parametric study of the barrier thickness (5-60 µm); the glassy
PET (gPET) diffusivity is given for reference.

The values reported in the article were computed with SFPPy 1.9.3.
Run from the SFPPy root:  python examples/TI_AG6062/case_C6_diffusion_barrier.py

@project: SFPPy - Safe Food Packaging in Python
@license: MIT
===============================================================================
"""

# =============================================================================
# SFPPy Path Bootstrap (allows running from examples/TI_AG6062/ without installation)
# =============================================================================
import sys
from pathlib import Path
_sfppy_root = Path(__file__).resolve().parents[2]
if str(_sfppy_root) not in sys.path:
    sys.path.insert(0, str(_sfppy_root))

import numpy as np
# NumPy 2.0 compat
if not hasattr(np, 'trapz'):
    np.trapz = np.trapezoid

def _f(x):
    """Convert numpy scalar/array to Python float."""
    return np.asarray(x).flat[0]

from patankar.loadpubchem import migrant
from patankar.geometry import Packaging3D
import patankar.food as food
import patankar.layer as polymer
from patankar.migration import senspatankar as solver
from patankar.migration import CFSimulationContainer as store

# %% 1) Geometry — 1 L cylinder
geom = Packaging3D('Cylinder', length=(19, "cm"), radius=(30, "mm"))
V_F, A = geom.get_volume_and_area()

# %% 2) Contact conditions
contactTemperature = (40, "degC")
contactTime = (10, "days")

# %% 3) Substance — toluene (small, fast-diffusing contaminant from recycling)
m = migrant("toluene")
print(f"Substance: toluene, MW = {_f(m.M):.0f} g/mol")
C0 = 200  # mg/kg in rPP layer
SML = 0.6  # mg/kg food (EU group SML for aromatic hydrocarbons)

# %% 4) Food — fatty
class fattyfood(food.realfood, food.semisolid, food.fat):
    name = "fatty food"

FOODlayer = fattyfood(
    volume=V_F, surfacearea=A,
    contacttime=contactTime, contacttemperature=contactTemperature,
    substance=m, simulant="ethanol"
)

# %% 5) WITHOUT barrier — monolayer rPP (300 µm)
rPP_mono = polymer.PP(l=(300, "um"), substance=m, C0=C0, T=contactTemperature)
sim_nofb = solver(rPP_mono, FOODlayer, name="rPP-no-FB")
CF_nofb = _f(sim_nofb.CFtarget)
s_nofb = 100 * CF_nofb / SML
print(f"\n=== WITHOUT functional barrier ===")
print(f"C_F = {CF_nofb:.3f} mg/kg food")
print(f"Severity = {s_nofb:.1f}% → {'PASS' if s_nofb<=100 else 'FAIL'}")

# %% 6) Barrier material choice — plasticized vs glassy PET
# IMPORTANT (two points):
#  (1) The FB layer must carry substance=m, otherwise SFPPy assigns a default
#      diffusivity (1e-14) instead of the hFV (hole-free-volume) prediction.
#  (2) The layer CLASS selects the right hFV parameterization automatically
#      (pure OO): gPET = glassy PET (Tg~76 C, dry virgin); wPET = wet/plasticized
#      PET (Tg~46 C). Recycled PET — or PET stored in a humid atmosphere — is
#      plasticized by sorbed water/contaminants, which depresses Tg and raises D.
#      We therefore model the realistic recyclate barrier with wPET, and quote
#      glassy gPET only as a (near-perfect) reference bound.
# hFV model: Zhu, Welle & Vitrac, Soft Matter 2019, 15(42), 8912-8932.

PET_glassy = polymer.gPET(l=(30, "um"), substance=m, C0=0, T=contactTemperature)
PET_plast  = polymer.wPET(l=(30, "um"), substance=m, C0=0, T=contactTemperature)
print(f"\n=== hFV diffusivity of toluene in PET at 40 C (auto-selected model) ===")
print(f"  gPET (glassy, Tg~{PET_glassy.Tg[0]:.0f} C, virgin ref): D = {_f(PET_glassy.D):.3e} m2/s")
print(f"  wPET (plasticized, Tg~{PET_plast.Tg[0]:.0f} C, recyclate): D = {_f(PET_plast.D):.3e} m2/s")
print(f"  ratio D(wPET)/D(gPET) = {_f(PET_plast.D)/_f(PET_glassy.D):.0f}x")

# %% 7) Parametric study — vary plasticized-PET (wPET) barrier thickness
#       Report lag time t_lag = l[0]^2 / (6 D[0]) of the food-side layer.
print(f"\n=== Parametric study: wPET (plasticized) barrier thickness, 10 d / 40 C ===")
print(f"{'wPET (µm)':>10} | {'t_lag (d)':>10} | {'t_lag/tc':>9} | {'C_F (mg/kg)':>12} | {'Severity':>10} | {'Reduction':>10} | {'Status':>8}")
print("-" * 90)
tc_s = 10 * 86400.0
thicknesses = [0, 5, 10, 15, 20, 30, 40, 50, 60]
for th in thicknesses:
    if th == 0:
        CF = CF_nofb
        tlag_d = None
    else:
        fb = polymer.wPET(l=(th, "um"), substance=m, C0=0, T=contactTemperature)
        rpp = polymer.PP(l=(300, "um"), substance=m, C0=C0, T=contactTemperature)
        bl = fb + rpp
        l0 = _f(np.asarray(bl.l).ravel()[0]); D0 = _f(np.asarray(bl.D).ravel()[0])
        tlag_d = l0**2 / (6 * D0) / 86400.0
        fl = fattyfood(volume=V_F, surfacearea=A,
                       contacttime=contactTime, contacttemperature=contactTemperature,
                       substance=m, simulant="ethanol")
        sim = solver(bl, fl, name=f"wPET{th}-rPP")
        CF = _f(sim.CFtarget)
    s = 100 * CF / SML
    red = CF_nofb / CF if CF > 0 else float('inf')
    status = "PASS" if s <= 100 else "FAIL"
    tlag_str = "-" if tlag_d is None else f"{tlag_d:.2f}"
    ratio_str = "-" if tlag_d is None else f"{tlag_d*86400.0/tc_s:.2f}"
    print(f"{th:>10} | {tlag_str:>10} | {ratio_str:>9} | {CF:>12.3f} | {s:>8.1f}% | {red:>8.1f}x | {status:>8}")

# %% 8) Key message
print(f"""
=== KEY MESSAGE ===
Toluene (MW=92, fast diffuser) in rPP (300 µm), 10 d / 40 C, SML = {SML} mg/kg:
- Without FB: C_F = {CF_nofb:.2f} mg/kg -> s = {s_nofb:.0f}% (FAIL).
- A plasticized-PET (wPET) functional barrier gives a GRADED reduction with
  thickness (no single magic value): ~20 um is the minimum for compliance here.
- The hFV diffusivity is set automatically by the layer class: glassy gPET
  (D~2.3e-18, near-perfect but unrealistic for a humid recyclate) vs wet wPET
  (D~4.9e-16). Recycled/humid PET is plasticized -> use wPET.
- Lag time t_lag = l[0]^2/(6 D[0]) of the barrier sets the deferral: thin layers
  (t_lag << t_contact) let toluene through; t_lag approaches t_contact near 60 um.
- Over long shelf life, even a compliant barrier eventually breaks through.
""")
