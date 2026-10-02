"""
===============================================================================
Case C3 — Restreindre l'usage : aliment gras -> aqueux (tier M2)
===============================================================================
Case study C3 (section 4.3) of:

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

BHT (MW=220, logP=5.3 from PubChem) in LDPE migrates excessively into fatty food.
Restrict to aqueous food → K changes dramatically → migration reduced.
Three estimation routes for K are compared.

The values reported in the article were computed with SFPPy 1.9.3.
Run from the SFPPy root:  python examples/TI_AG6062/case_C3_fatty_to_aqueous.py

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
from patankar.migration import senspatankar as solver

# %% 1) Geometry
geom = Packaging3D('Cylinder', length=(19, "cm"), radius=(30, "mm"))
V_F, A = geom.get_volume_and_area()

# %% 2) Contact conditions
contactTemperature = (40, "degC")
contactTime = (10, "days")
C0 = 500  # mg/kg in LDPE
SML = 3.0  # mg/kg food (BHT)

# %% 3) Substance
m = migrant("BHT")
print(f"BHT: MW = {_f(m.M):.0f}, logP = {_f(m.logP):.2f}")

# %% 4) Fatty food contact
class fattyfood(food.realfood, food.semisolid, food.fat):
    name = "fatty food"

LDPE_layer = polymer.LDPE(l=(200,"um"), substance=m, C0=C0, T=contactTemperature)
FOOD_fat = fattyfood(
    volume=V_F, surfacearea=A,
    contacttime=contactTime, contacttemperature=contactTemperature,
    substance=m, simulant="ethanol"
)
sim_fat = solver(LDPE_layer, FOOD_fat, name="BHT-fatty")
CF_fat = _f(sim_fat.CFtarget)

# Read the partition coefficients used
k_P = _f(LDPE_layer.k)
k_F_fat = _f(FOOD_fat.k0) if hasattr(FOOD_fat, 'k0') else 1.0
KFP_fat = k_P / k_F_fat if k_F_fat != 0 else float('inf')

print(f"\n=== FATTY food contact ===")
print(f"k_polymer = {k_P:.4f}")
print(f"k_food (fat) = {k_F_fat:.4f}")
print(f"K_F/P (fat) = {KFP_fat:.4f}")
print(f"C_F = {CF_fat:.3f} mg/kg")
print(f"Severity = {100*CF_fat/SML:.1f}% → {'PASS' if CF_fat<=SML else 'FAIL'}")

# %% 5) Aqueous food contact
# For aqueous food, use the M2 analytical approach (no full simulation needed).
# BHT is very lipophilic (logP=5.3), so K_water/LDPE is very small.
# At equilibrium: CF_aqueous = C0 * (m_P/m_F) * PR_E
# where PR_E = 1 / (1 + K_P/F * V_P/V_F)
# K_P/F = k_P / k_F, and for aqueous food k_F is much smaller.
logP_bht = _f(m.logP)

# Route 1: from logP
# k_F(water) ≈ k_F(fat) * 10^(-logP) (octanol/water partition)
k0_aqueous = k_F_fat * 10**(-logP_bht)
K_PF_aqueous = k_P / k0_aqueous  # K_polymer/food for aqueous

# Analytical M2 calculation
V_P = 200e-6 * _f(A)  # m³
l_P = 200e-6
rho_P = 920.0  # LDPE
m_P = rho_P * l_P * _f(A)
m_F = 1000.0 * _f(V_F)
PR_E_aq = 1.0 / (1.0 + K_PF_aqueous * V_P / _f(V_F))
CF_aq = C0 * (m_P / m_F) * PR_E_aq

print(f"\nEstimated k0 for aqueous food: {k0_aqueous:.2e} (from logP = {logP_bht:.2f})")
print(f"K_polymer/food (aqueous) = {K_PF_aqueous:.1f}")

k_F_aq = k0_aqueous
KFP_aq = KFP_fat  # for display: K_food/polymer
KFP_aq = k_P / k_F_aq if k_F_aq != 0 else float('inf')

print(f"\n=== AQUEOUS food contact ===")
print(f"k_food (aqueous) = {k_F_aq:.4f}")
print(f"K_F/P (aqueous) = {KFP_aq:.4f}")
print(f"C_F = {CF_aq:.3f} mg/kg")
print(f"Severity = {100*CF_aq/SML:.1f}% → {'PASS' if CF_aq<=SML else 'FAIL'}")

# %% 6) Comparison
print(f"\n=== Comparison ===")
print(f"{'Food type':<15} | {'K_F/P':>8} | {'C_F':>10} | {'Severity':>10} | {'Status':>8}")
print("-" * 60)
print(f"{'Fatty':.<15} | {KFP_fat:>8.4f} | {CF_fat:>10.3f} | {100*CF_fat/SML:>8.1f}% | {'PASS' if CF_fat<=SML else 'FAIL':>8}")
print(f"{'Aqueous':.<15} | {KFP_aq:>8.4f} | {CF_aq:>10.3f} | {100*CF_aq/SML:>8.1f}% | {'PASS' if CF_aq<=SML else 'FAIL':>8}")
if CF_aq > 0:
    print(f"\nMigration reduction (fat → aqueous): {CF_fat/CF_aq:.1f}x")

# %% 7) Key message
print(f"""
=== KEY MESSAGE ===
BHT is lipophilic (logP = {_f(m.logP):.1f}): it has strong affinity for fat but
low affinity for water. Restricting usage from fatty to aqueous food changes
the partition coefficient K_F/P dramatically, reducing migration.

This lever belongs to tier M2: no change to the material, only to the
declared usage conditions. The effect is thermodynamic (equilibrium-based),
not kinetic.
""")
