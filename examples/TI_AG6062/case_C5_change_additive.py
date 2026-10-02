"""
===============================================================================
Case C5 — Changer d'additif (tier M3, recalcul de D)
===============================================================================
Case study C5 (section 4.5) of:

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

Replace a small, mobile antioxidant (BHT, MW=220) with a larger, slower one
(Irganox 1076, MW=531). Both are phenolic antioxidants. The Piringer model
predicts D ~ exp(-0.1351 * M^(2/3)), so the larger molecule diffuses much
more slowly.

The values reported in the article were computed with SFPPy 1.9.3.
Run from the SFPPy root:  python examples/TI_AG6062/case_C5_change_additive.py

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
# NumPy 2.0 compat: np.trapz was removed, aliased to np.trapezoid
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

# %% 1) Geometry — 1 L cylindrical container
geom = Packaging3D('Cylinder', length=(19, "cm"), radius=(30, "mm"))
V_F, A = geom.get_volume_and_area()

# %% 2) Contact conditions
contactTemperature = (40, "degC")  # EU regulatory standard
contactTime = (10, "days")

# %% 3) Define both additives
m_bht = migrant("BHT")           # MW ~ 220
m_irg = migrant("Irganox 1076")  # MW ~ 531
print(f"BHT:          MW = {_f(m_bht.M):.0f} g/mol")
print(f"Irganox 1076: MW = {_f(m_irg.M):.0f} g/mol")

# %% 4) Define LDPE layers
C0 = 500  # mg/kg (same mass concentration for both)
LDPE_bht = polymer.LDPE(
    l=(200, "um"),
    substance=m_bht,
    C0=C0,
    T=contactTemperature
)
LDPE_irg = polymer.LDPE(
    l=(200, "um"),
    substance=m_irg,
    C0=C0,
    T=contactTemperature
)

# Print diffusivities from Piringer model
print(f"\n=== Diffusivities (Piringer model, LDPE @ 40°C) ===")
print(f"D(BHT)          = {_f(LDPE_bht.D):.3e} m²/s")
print(f"D(Irganox 1076) = {_f(LDPE_irg.D):.3e} m²/s")
print(f"Ratio D(BHT)/D(I1076) = {_f(LDPE_bht.D)/_f(LDPE_irg.D):.1f}x")

# %% 5) Define food — fatty food
class fattyfood(food.realfood, food.semisolid, food.fat):
    name = "fatty food"

FOODlayer = fattyfood(
    volume=V_F,
    surfacearea=A,
    contacttime=contactTime,
    contacttemperature=contactTemperature,
    substance=m_bht,
    simulant="ethanol"
)

# %% 6) Run migration — BHT
sim_bht = solver(LDPE_bht, FOODlayer, name="BHT-LDPE")
CF_bht = _f(sim_bht.CFtarget)
print(f"\n=== Migration results (10 d @ 40°C) ===")
print(f"BHT:          C_F = {CF_bht:.3f} mg/kg food")

# %% 7) Run migration — Irganox 1076
FOODlayer_irg = fattyfood(
    volume=V_F,
    surfacearea=A,
    contacttime=contactTime,
    contacttemperature=contactTemperature,
    substance=m_irg,
    simulant="ethanol"
)
sim_irg = solver(LDPE_irg, FOODlayer_irg, name="I1076-LDPE")
CF_irg = _f(sim_irg.CFtarget)
print(f"Irganox 1076: C_F = {CF_irg:.3f} mg/kg food")

# %% 8) Severity comparison
SML_bht = 3.0  # mg/kg food (EU 10/2011)
SML_irg = 6.0  # mg/kg food (EU 10/2011)
s_bht = 100 * CF_bht / SML_bht
s_irg = 100 * CF_irg / SML_irg

print(f"\n=== Severity ===")
print(f"{'Substance':<15} | {'MW':>5} | {'D (m²/s)':>12} | {'C_F':>8} | {'SML':>6} | {'s (%)':>8} | {'Status':>8}")
print("-" * 80)
print(f"{'BHT':<15} | {_f(m_bht.M):>5.0f} | {_f(LDPE_bht.D):>12.3e} | {CF_bht:>8.3f} | {SML_bht:>6.1f} | {s_bht:>8.1f} | {'PASS' if s_bht<=100 else 'FAIL':>8}")
print(f"{'Irganox 1076':<15} | {_f(m_irg.M):>5.0f} | {_f(LDPE_irg.D):>12.3e} | {CF_irg:>8.3f} | {SML_irg:>6.1f} | {s_irg:>8.1f} | {'PASS' if s_irg<=100 else 'FAIL':>8}")
print(f"\nMigration reduction factor: {CF_bht/CF_irg:.1f}x")

# %% 9) Fourier numbers — explain why migration is similar
D_bht = _f(LDPE_bht.D)
D_irg = _f(LDPE_irg.D)
l = 200e-6  # m
t = 10 * 86400  # s
Fo_bht = D_bht * t / l**2
Fo_irg = D_irg * t / l**2
print(f"\n=== Fourier numbers (Fo = D*t/l²) ===")
print(f"Fo(BHT)     = {Fo_bht:.2f}  {'(>> 1: near equilibrium)' if Fo_bht > 1 else '(< 1: kinetic regime)'}")
print(f"Fo(I1076)   = {Fo_irg:.2f}  {'(>> 1: near equilibrium)' if Fo_irg > 1 else '(< 1: kinetic regime)'}")

# %% 10) Shorter contact — 1 day @ 25°C to reveal kinetic difference
print(f"\n=== Shorter contact: 1 day @ 25°C ===")
T_short = (25, "degC")
t_short = (1, "days")
LDPE_bht_short = polymer.LDPE(l=(200,"um"), substance=m_bht, C0=C0, T=T_short)
LDPE_irg_short = polymer.LDPE(l=(200,"um"), substance=m_irg, C0=C0, T=T_short)
FOOD_short = fattyfood(volume=V_F, surfacearea=A, contacttime=t_short,
                        contacttemperature=T_short, substance=m_bht, simulant="ethanol")
FOOD_short_irg = fattyfood(volume=V_F, surfacearea=A, contacttime=t_short,
                            contacttemperature=T_short, substance=m_irg, simulant="ethanol")
sim_bht_short = solver(LDPE_bht_short, FOOD_short, name="BHT-1d-25C")
sim_irg_short = solver(LDPE_irg_short, FOOD_short_irg, name="I1076-1d-25C")
CF_bht_s = _f(sim_bht_short.CFtarget)
CF_irg_s = _f(sim_irg_short.CFtarget)
print(f"BHT:          C_F(1d@25C) = {CF_bht_s:.3f} mg/kg")
print(f"Irganox 1076: C_F(1d@25C) = {CF_irg_s:.3f} mg/kg")
if CF_irg_s > 0:
    print(f"Reduction factor: {CF_bht_s/CF_irg_s:.1f}x")

# %% 11) Key message
print(f"""
=== KEY MESSAGE ===
Replacing BHT (MW=220) by Irganox 1076 (MW=531) in LDPE:
- Diffusivity decreases by {D_bht/D_irg:.0f}x (Piringer model)
- At 10d @ 40°C: both near equilibrium (Fo >> 1) → similar migration
- At 1d @ 25°C: kinetic regime → migration reduced by {CF_bht_s/CF_irg_s:.0f}x
- The kinetic advantage of larger molecules is most effective for:
  (a) thicker films, (b) shorter contacts, (c) lower temperatures
- For thin films at long contact, the PARTITION COEFFICIENT (M2) matters
  more than diffusivity (M3)
""")
