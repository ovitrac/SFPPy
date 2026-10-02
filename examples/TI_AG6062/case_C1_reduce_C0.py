"""
===============================================================================
Case C1 — Diminuer C0 (tier M0)
===============================================================================
Case study C1 (section 4.1) of:

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

Demonstrates the simplest design lever: reducing initial concentration C0.
A substance (Irganox 1076) in recycled PP at C0=1000 mg/kg exceeds the SML
at tier M0 (total transfer). We show:
  1. The M0 mass balance (total transfer → severity)
  2. The threshold C0_max for compliance
  3. The effect of blending virgin + recycled PP

The values reported in the article were computed with SFPPy 1.9.3.
Run from the SFPPy root:  python examples/TI_AG6062/case_C1_reduce_C0.py

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

# %% 1) Define geometry — 1 L cylindrical container
geom = Packaging3D('Cylinder', length=(19, "cm"), radius=(30, "mm"))
V_F, A = geom.get_volume_and_area()  # m³, m²
V_F, A = np.asarray(V_F).item(), np.asarray(A).item()
print(f"Geometry: V_F = {V_F*1e6:.1f} mL, A = {A*1e4:.1f} cm²")

# %% 2) Define substance
m = migrant("Irganox 1076")
print(f"Substance: {m.name[6]}, MW = {float(m.M):.0f} g/mol")

# SML for Irganox 1076: 6 mg/kg food (EU 10/2011)
SML = 6.0  # mg/kg food

# %% 3) Material parameters (rPP monolayer)
l_P = 300e-6        # 300 µm thickness [m]
rho_P = 900.0       # PP density [kg/m³]
rho_F = 1000.0      # food density [kg/m³]
C0_initial = 1000.0 # mg/kg polymer

# Mass of polymer and food
m_P = rho_P * l_P * A    # kg
m_F = rho_F * V_F         # kg
print(f"Masses: m_P = {m_P*1e3:.2f} g, m_F = {m_F*1e3:.1f} g")

# %% 4) Tier M0: total transfer
# All substance migrates into food
C_F_M0 = C0_initial * m_P / m_F  # mg/kg food
s_M0 = 100 * C_F_M0 / SML
print(f"\n=== Tier M0 (total transfer) ===")
print(f"C0 = {C0_initial} mg/kg polymer")
print(f"C_F(M0) = {C_F_M0:.2f} mg/kg food")
print(f"Severity s = {s_M0:.1f} %")
print(f"Compliant: {'YES' if s_M0 <= 100 else 'NO'}")

# %% 5) Threshold C0_max
C0_max = SML * m_F / m_P  # mg/kg polymer
print(f"\nC0_max for compliance at M0 = {C0_max:.1f} mg/kg polymer")

# %% 6) Blending virgin + recycled
print(f"\n=== Blending strategy (M0) ===")
print(f"{'Recycled %':>12} | {'C0_eff':>10} | {'C_F(M0)':>10} | {'Severity':>10} | {'Status':>8}")
print("-" * 65)
C0_recycled = 1000.0
for x_recycled in [1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]:
    C0_eff = x_recycled * C0_recycled
    CF = C0_eff * m_P / m_F
    s = 100 * CF / SML
    status = "PASS" if s <= 100 else "FAIL"
    print(f"{x_recycled*100:>10.0f} % | {C0_eff:>8.0f}   | {CF:>8.2f}   | {s:>8.1f} % | {status:>8}")

# %% 7) Key message
print(f"""
=== KEY MESSAGE ===
At tier M0 (total transfer, most conservative):
- C_F is proportional to C0 × (m_P/m_F)
- Ratio m_P/m_F = {m_P/m_F:.4f}
- C0_max = {C0_max:.1f} mg/kg for SML = {SML} mg/kg
- Even at 30% recycled content, compliance may require C0 < {C0_max:.0f} mg/kg
- M0 is ultra-conservative; higher tiers (M1, M2, M3) will show lower migration
""")
