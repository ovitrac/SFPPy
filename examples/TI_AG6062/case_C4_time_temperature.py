"""
===============================================================================
Case C4 — Modifier le parcours temps/température (tier M3)
===============================================================================
Case study C4 (section 4.4) of:

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

Irganox 1076 in PP (300 µm, C0 = 500 mg/kg) in contact with a fatty food.
Regulatory test: 10 d at 40 °C (compliant); real use: 2 d at 25 °C (chilled
ready-meal). Shows the Arrhenius effect D(T): the SML (6 mg/kg) is exceeded
from 50 °C. Section 6b reproduces table 16 of the article (2, 10, 30 d).

The values reported in the article were computed with SFPPy 1.9.3.
Run from the SFPPy root:  python examples/TI_AG6062/case_C4_time_temperature.py

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

# %% 2) Substance
m = migrant("Irganox 1076")
C0 = 500
SML = 6.0  # mg/kg food

# %% 3) Regulatory conditions: 10d @ 40C
T_reg, t_reg = (40, "degC"), (10, "days")

class fattyfood(food.realfood, food.semisolid, food.fat):
    name = "fatty food"

PP_reg = polymer.PP(l=(300, "um"), substance=m, C0=C0, T=T_reg)
FOOD_reg = fattyfood(volume=V_F, surfacearea=A,
                     contacttime=t_reg, contacttemperature=T_reg,
                     substance=m, simulant="ethanol")
sim_reg = solver(PP_reg, FOOD_reg, name="I1076-PP-40C-10d")
CF_reg = _f(sim_reg.CFtarget)
D_reg = _f(PP_reg.D)

print(f"=== Regulatory: 10 d @ 40°C ===")
print(f"D(40°C) = {D_reg:.3e} m²/s")
print(f"C_F = {CF_reg:.3f} mg/kg → s = {100*CF_reg/SML:.1f}% → {'PASS' if CF_reg<=SML else 'FAIL'}")

# %% 4) Real use: 2d @ 25C
T_real, t_real = (25, "degC"), (2, "days")

PP_real = polymer.PP(l=(300, "um"), substance=m, C0=C0, T=T_real)
FOOD_real = fattyfood(volume=V_F, surfacearea=A,
                      contacttime=t_real, contacttemperature=T_real,
                      substance=m, simulant="ethanol")
sim_real = solver(PP_real, FOOD_real, name="I1076-PP-25C-2d")
CF_real = _f(sim_real.CFtarget)
D_real = _f(PP_real.D)

print(f"\n=== Real use: 2 d @ 25°C ===")
print(f"D(25°C) = {D_real:.3e} m²/s")
print(f"C_F = {CF_real:.3f} mg/kg → s = {100*CF_real/SML:.1f}% → {'PASS' if CF_real<=SML else 'FAIL'}")

# %% 5) Fourier numbers
l = 300e-6
Fo_reg = D_reg * 10*86400 / l**2
Fo_real = D_real * 2*86400 / l**2
print(f"\n=== Fourier numbers ===")
print(f"Fo(40°C, 10d) = {Fo_reg:.4f}")
print(f"Fo(25°C, 2d)  = {Fo_real:.4f}")

# %% 6) Multi-temperature comparison
print(f"\n=== Parametric: same 10 d, varying temperature ===")
print(f"{'T (°C)':>8} | {'D (m²/s)':>12} | {'Fo':>8} | {'C_F':>10} | {'s (%)':>8} | {'Status':>8}")
print("-" * 70)
for T_val in [5, 10, 15, 20, 25, 30, 35, 40, 50, 60]:
    T_i = (T_val, "degC")
    PP_i = polymer.PP(l=(300, "um"), substance=m, C0=C0, T=T_i)
    F_i = fattyfood(volume=V_F, surfacearea=A,
                    contacttime=(10, "days"), contacttemperature=T_i,
                    substance=m, simulant="ethanol")
    sim_i = solver(PP_i, F_i, name=f"T{T_val}")
    CF_i = _f(sim_i.CFtarget)
    D_i = _f(PP_i.D)
    Fo_i = D_i * 10*86400 / l**2
    s_i = 100 * CF_i / SML
    print(f"{T_val:>8} | {D_i:>12.3e} | {Fo_i:>8.4f} | {CF_i:>10.3f} | {s_i:>8.1f} | {'PASS' if s_i<=100 else 'FAIL':>8}")

# %% 6b) Table 16 of the article: temperature x contact time (2, 10, 30 d)
print(f"\n=== Table 16: C_F (mg/kg) versus temperature and contact time ===")
print(f"{'T (°C)':>8} | {'D (m²/s)':>12} | {'Fo (10 d)':>9} | {'2 d':>7} | {'10 d':>7} | {'30 d':>7}")
print("-" * 66)
for T_val in [5, 25, 40, 50, 60]:
    row = []
    for t_days in (2, 10, 30):
        PP_i = polymer.PP(l=(300, "um"), substance=m, C0=C0, T=(T_val, "degC"))
        F_i = fattyfood(volume=V_F, surfacearea=A,
                        contacttime=(t_days, "days"), contacttemperature=(T_val, "degC"),
                        substance=m, simulant="ethanol")
        CF_i = _f(solver(PP_i, F_i, name=f"T{T_val}-{t_days}d").CFtarget)
        row.append(f"{CF_i:.3g}{'' if CF_i <= SML else '*'}")
    D_i = _f(PP_i.D)
    print(f"{T_val:>8} | {D_i:>12.3e} | {D_i*10*86400/l**2:>9.4f} | " + " | ".join(f"{c:>7}" for c in row))
print("(*) non-compliant: C_F > SML = 6 mg/kg")

# %% 7) Key message
print(f"""
=== KEY MESSAGE ===
D follows Arrhenius: D(T) = D0 * exp(-Ea/RT)
- D(40°C) / D(25°C) = {D_reg/D_real:.1f}x
- Real-use conditions (2d@25°C) vs regulatory (10d@40°C):
  Migration reduced from {CF_reg:.2f} to {CF_real:.3f} mg/kg ({CF_reg/max(CF_real,1e-10):.0f}x)
- Temperature is a powerful lever: lowering T from 40 to 25°C reduces
  D by {D_reg/D_real:.0f}x, and combined with shorter contact time, migration
  drops dramatically
""")
