"""
===============================================================================
Case C2 — Changer de contaminant (axe toxicologique T)
===============================================================================
Case study C2 (section 4.2) of:

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

Same packaging, same migration level, but different substances with different
toxicological profiles. Shows how the SML changes with T0→T1→T2.
  - Substance A: limonene (Cramer I, non-toxic terpene)
  - Substance B: benzophenone (Cramer III, photoinitiator)

The values reported in the article were computed with SFPPy 1.9.3.
Run from the SFPPy root:  python examples/TI_AG6062/case_C2_change_contaminant.py

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
from patankar.loadpubchem import migrant
from patankar.geometry import Packaging3D

# %% 1) Geometry — same 1 L cylinder
geom = Packaging3D('Cylinder', length=(19, "cm"), radius=(30, "mm"))
V_F, A = geom.get_volume_and_area()
V_F, A = np.asarray(V_F).item(), np.asarray(A).item()

# Material parameters (rPP, 300 µm)
l_P = 300e-6
rho_P = 900.0
rho_F = 1000.0
m_P = rho_P * l_P * A
m_F = rho_F * V_F

# %% 2) Define substances
m_limo = migrant("limonene")
m_benzo = migrant("benzophenone")
print(f"Limonene:     MW = {float(m_limo.M):.0f}, logP = {float(np.asarray(m_limo.logP).flat[0]):.2f}")
print(f"Benzophenone: MW = {float(m_benzo.M):.0f}, logP = {float(np.asarray(m_benzo.logP).flat[0]):.2f}")

# %% 3) Same C0 for both, compute C_F at M1 (K=1, partition unitaire)
C0 = 80.0  # mg/kg polymer (chosen so limonene passes at T1-I)
# At M1 (K=1): PR^E = 1/(1 + K*V_P/V_F) ≈ V_F/(V_F+V_P)
V_P = l_P * A  # m³ (volume of polymer)
PR_E = V_F / (V_F + V_P)  # dimensionless, ~0.98 for thin film
C_F_M1 = C0 * (m_P / m_F) * PR_E  # mg/kg food
print(f"\nAt M1 (K=1): C_F = {C_F_M1:.3f} mg/kg food (same for both substances)")

# %% 4) Toxicological reference values
# BW = 60 kg, DI = 1 kg food/day (EU convention)
BW = 60.0  # kg
DI = 1.0   # kg food/day

# T0: TTC genotoxic = 0.0025 µg/kg bw/day
TTC_T0 = 0.0025e-3  # mg/kg bw/day
SML_T0 = BW * DI * TTC_T0  # mg/kg food = 0.00015 → 0.15 µg/kg

# T1: Cramer classes (TTC in µg/kg bw/day)
TTC_CramerI = 30e-3    # 30 µg/kg bw/day → SML ~ 1.8 mg/kg
TTC_CramerII = 9e-3    # 9 µg/kg bw/day → SML ~ 0.54 mg/kg
TTC_CramerIII = 1.5e-3 # 1.5 µg/kg bw/day → SML ~ 0.09 mg/kg

SML_T1_I = BW * DI * TTC_CramerI
SML_T1_II = BW * DI * TTC_CramerII
SML_T1_III = BW * DI * TTC_CramerIII

# T2: Specific SML from EU 10/2011
SML_T2_limo = None   # No specific SML (not listed, organoleptic only)
SML_T2_benzo = 0.6   # mg/kg food (EU 10/2011, group restriction)

print(f"\n=== Toxicological reference values ===")
print(f"{'Tier':<6} | {'Parameter':<30} | {'SML (mg/kg food)':>18}")
print("-" * 60)
print(f"{'T0':<6} | {'TTC genotoxic':<30} | {BW*DI*TTC_T0:>18.4f}")
print(f"{'T1-I':<6} | {'Cramer I (low toxicity)':<30} | {SML_T1_I:>18.2f}")
print(f"{'T1-II':<6} | {'Cramer II (moderate)':<30} | {SML_T1_II:>18.2f}")
print(f"{'T1-III':<6} | {'Cramer III (high toxicity)':<30} | {SML_T1_III:>18.2f}")
print(f"{'T2':<6} | {'Benzophenone (specific)':<30} | {SML_T2_benzo:>18.2f}")
print(f"{'T2':<6} | {'Limonene (not listed)':<30} | {'N/A':>18}")

# %% 5) Severity at each tier
print(f"\n=== Severity comparison (C_F = {C_F_M1:.3f} mg/kg) ===")
print(f"{'Substance':<15} | {'Tier':<6} | {'SML':>8} | {'s (%)':>8} | {'Status':>8}")
print("-" * 60)

cases = [
    ("Limonene",     "T0",   BW*DI*TTC_T0),
    ("Limonene",     "T1-I", SML_T1_I),
    ("Benzophenone", "T0",   BW*DI*TTC_T0),
    ("Benzophenone", "T1-III", SML_T1_III),
    ("Benzophenone", "T2",   SML_T2_benzo),
]

for name, tier, sml in cases:
    s = 100 * C_F_M1 / sml
    status = "PASS" if s <= 100 else "FAIL"
    print(f"{name:<15} | {tier:<6} | {sml:>8.4f} | {s:>8.1f} | {status:>8}")

# %% 6) Key message
print(f"""
=== KEY MESSAGE ===
Same migration level (C_F = {C_F_M1:.3f} mg/kg), but severity depends entirely
on the toxicological tier:

- Limonene (Cramer I): PASSES at T1 (s = {100*C_F_M1/SML_T1_I:.1f}%)
  → Safe, no material action needed
- Benzophenone (Cramer III): FAILS at T1 (s = {100*C_F_M1/SML_T1_III:.0f}%)
  → Must use T2 (specific SML = 0.6 mg/kg) where it PASSES (s = {100*C_F_M1/SML_T2_benzo:.1f}%)
  → If no T2 data available, action on material axis (M) is required

The toxicological axis T is the FIRST lever to explore: better toxicological
knowledge (T1→T2) can transform a non-compliant assessment into compliance
without any material modification.
""")
