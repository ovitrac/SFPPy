"""
===============================================================================
Case C7 — Barrière fonctionnelle : diffusivité et solubilité couplées (tier M2+M3 + FB)
===============================================================================
Case study C7 (section 4.7) of:

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

Limonene (non-polar, logP = 3.4 from PubChem) migrates from a recycled PP (rPP, 300 um,
500 mg/kg) into a fatty food. Three virgin functional barriers (20 um, C0 = 0)
are compared: virgin PP (same material: the only "pure D" barrier, acting
through thickness alone), PA6 (polar, but not polar enough to exclude
limonene) and PET (glassy, diffusion barrier). The PA6 effect is split into
its solubility part (K only) and its diffusion part (D only) with layerLink.

Physics: a thin barrier acts through its permeance D_j K_{j/s} / l_j. It delays
migration but cannot change the equilibrium, which is set by the mass balance
of the source and the food: every barrier effect is transient.

The values reported in the article were computed with SFPPy 1.9.3.
Run from the SFPPy root:  python examples/TI_AG6062/case_C7_functional_barrier.py

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
if not hasattr(np, 'trapz'):
    np.trapz = np.trapezoid

def _f(x):
    return np.asarray(x).flat[0]

from patankar.loadpubchem import migrant
from patankar.geometry import Packaging3D
import patankar.food as food
import patankar.layer as polymer
from patankar.layer import layerLink
from patankar.migration import senspatankar as solver

# %% 1) Geometry and conditions
geom = Packaging3D('Cylinder', length=(19, "cm"), radius=(30, "mm"))
V_F, A = geom.get_volume_and_area()
T = (40, "degC")
TIMES = [10, 100, 1000, 10000]          # days
m = migrant("limonene")
C0 = 500                                 # mg/kg in rPP
SML = 1.0                                # mg/kg food (organoleptic threshold)
print(f"Limonene: MW = {_f(m.M):.0f}, logP = {_f(m.logP):.2f}")

class fattyfood(food.realfood, food.semisolid, food.fat):
    name = "fatty food"

def food_at(days):
    return fattyfood(volume=V_F, surfacearea=A, contacttime=(days, "days"),
                     contacttemperature=T, substance=m, simulant="ethanol")

def source():
    return polymer.PP(l=(300, "um"), substance=m, C0=C0, T=T)

def barrier(cls):
    return cls(l=(20, "um"), substance=m, C0=0, T=T)   # virgin barrier

# %% 2) Henry-like coefficients and crystallinity
PP_ref, PA6_ref = source(), barrier(polymer.PA6)
k_PP, k_PA6 = _f(PP_ref.k), _f(PA6_ref.k)
c_PP, c_PA6 = PP_ref.crystallinity_history[0], PA6_ref.crystallinity_history[0]
crys = (1 - c_PP) / (1 - c_PA6)
print(f"\nk_PA6 / k_PP = {k_PA6 / k_PP:.2f} = {k_PA6 / k_PP / crys:.2f} (amorphous phases)"
      f" x (1 - c_PP)/(1 - c_PA6) = {crys:.3f}   [c_PP = {c_PP}, c_PA6 = {c_PA6}]")
D_PP, D_PA6 = _f(PP_ref.D), _f(PA6_ref.D)
print(f"D_PA6 / D_PP = {D_PA6 / D_PP:.3e}   (D_PP = {D_PP:.3e}, D_PA6 = {D_PA6:.3e} m2/s)")

# %% 3) Configurations
def red_str(r, width=8):
    """Reduction factor, capped for display (C_F numerically zero behind PET)."""
    return f"{'>1e6x':>{width + 1}}" if r > 1e6 else f"{r:>{width}.3g}x"

def config(name):
    if name == "rPP monolayer":
        return source()
    if name == "vPP(20um) + rPP":
        return barrier(polymer.PP) + source()
    if name == "PA6(20um) + rPP":
        return barrier(polymer.PA6) + source()
    if name == "PA6, K only":                  # D of the barrier forced to D_PP
        ml = barrier(polymer.PA6) + source()
        ml.Dlink = layerLink("D", indices=[0], values=[D_PP])
        return ml
    if name == "PA6, D only":                  # k of the barrier forced to k_PP
        ml = barrier(polymer.PA6) + source()
        ml.klink = layerLink("k", indices=[0], values=[k_PP])
        return ml
    if name == "PET(20um) + rPP":
        return barrier(polymer.gPET) + source()
    raise ValueError(name)

NAMES = ["rPP monolayer", "vPP(20um) + rPP", "PA6(20um) + rPP", "PA6, K only", "PA6, D only", "PET(20um) + rPP"]
CF = {n: {} for n in NAMES}
for n in NAMES:
    for t in TIMES:
        CF[n][t] = _f(solver(config(n), food_at(t), name=f"{n}-{t}d").CFtarget)

# %% 4) Results at 10 days
print(f"\n=== 10 days, 40 degC ===")
print(f"{'Configuration':<20} | {'C_F (mg/kg)':>11} | {'s (%)':>7} | {'Reduction':>9} | {'Status':>6}")
print("-" * 66)
for n in NAMES:
    cf = CF[n][10]; red = CF["rPP monolayer"][10] / cf if cf > 0 else float("inf")
    print(f"{n:<20} | {cf:>11.3g} | {100 * cf / SML:>7.1f} | {red_str(red)} | {'PASS' if cf <= SML else 'FAIL':>6}")

# %% 5) Transient nature of every barrier effect
print(f"\n=== Reduction factor C_F(rPP) / C_F(config) versus contact time ===")
print(f"{'Configuration':<20} | " + " | ".join(f"{t:>7d} d" for t in TIMES))
print("-" * 66)
for n in NAMES[1:]:
    print(f"{n:<20} | " + " | ".join(red_str(CF['rPP monolayer'][t] / CF[n][t]) if CF[n][t] > 0 else f"{'inf':>9}" for t in TIMES))
print(f"{'C_F rPP (mg/kg)':<20} | " + " | ".join(f"{CF['rPP monolayer'][t]:>9.3g}" for t in TIMES))

# %% 6) Key message
print("""
=== KEY MESSAGE ===
- No real barrier acts through D alone or K alone. A pure-D barrier exists
  only when the barrier is the same material as the source (virgin PP over
  rPP): it acts through thickness alone. A pure-K barrier does not exist:
  a more polar polymer is also more cohesive, more crystalline, with a higher
  Tg and therefore a lower D.
- PA6 is not polar enough to exclude limonene: k_PA6/k_PP stays close to 1
  once crystallinities are accounted for; its effect comes from D.
- A thin barrier acts through its permeance D_j K_{j/s} / l_j. It delays
  migration but cannot change the equilibrium set by the source and the food:
  every barrier effect, of solubility or of diffusion, fades at long times.
""")
