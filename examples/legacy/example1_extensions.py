#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example 1 — extensions
Simulation of additive transfer through a single LDPE layer (monolayer).

Objective
---------
Estimate the partition and potential migration of a
phenolic antioxidant (Irganox 1076) from a thin LDPE
film into a model sandwich food matrix, then evaluate
the effect of contact time and temperature using
simplified M2 (equilibrium) and M3 (kinetics) tiers.

Description
-----------
The example proceeds as follows:

1. Geometry:
   - Define a cylindrical packaging geometry
     (height = 19 cm, radius = 30 mm).
   - Compute internal volume and food contact area.

2. Migrant and Film:
   - Load Irganox 1076 properties from public
     databases via `migrant()`.
   - Define a 100 µm LDPE layer containing the migrant.

3. Food Proxy:
   - Define a semi-solid, fatty "sandwich" simulant
     (ethanol 95%) with specified contact time
     (10 days at 7 °C).

4. Equilibrium (Tier M2):
   - Compute the partition coefficient K_(F/LDPE).
   - Values < 1 indicate stronger affinity for the polymer.

5. Kinetics (Tier M3):
   - Perform migration calculation at 10 d / 7 °C.
   - Output concentration in food (mg/kg) and
     severity relative to the SML.

6. Additional Step:
   - Add 4 h at 25 °C (warming before consumption).
   - Chain simulations (`>>` operator).
   - Compute updated food concentration and severity.

Outputs
-------
- Internal geometry (volume, area).
- Migrant identity and properties.
- K_(F/LDPE) partition coefficient.
- Concentration in food (mg/kg).
- Severity relative to SML (%) at each step.
- Optional plots:
  - `plotCF()` : time-course in food
  - `plotCx()` : concentration profiles across layers

References
----------
- Tiered approach to migration modelling (M2:
  equilibrium; M3: diffusion-kinetics).
- Public databases for migrant properties (PubChem).

Notes
-----
This example demonstrates how SFPPy can carry out
basic Tier M2 and M3 assessments from minimal
inputs (compound name, geometry, film thickness).

Outputs
--------
[Geometry] V = 0.000537212 m^3
[Geometry] A = 0.041469 m^2
[Migrant] irganox 1076 | CAS=['2082-79-3'] | M=530.9 g/mol | logP=[13.8]
[M2] K_(F/LDPE) = 6.741 (<1: stronger LDPE affinity)
[M3] 10 d @ 7 C: C_F = 5.22691 mg/kg, severity = 87.12
[M3] +4 h @ 25 C: C_F = 5.56426 mg/kg, severity = 92.74

Created on Mon Aug 25 14:46:59 2025

@author: olivier.vitrac@gmail.com
Generative Simulation Initiative

"""

# --- Imports -------------------------------------------------------------
from patankar.layer import LDPE          # LDPE model
from patankar.loadpubchem import migrant # PubChem + props
from patankar.geometry import Packaging3D# 3D geometry
import patankar.food as food             # food/simulants

# --- 1) Geometry: cylinder (height x radius) -----------------------------
# Units as tuples (value, "unit"), handled by SFPPy.
internal_volume, surface_area = Packaging3D(
   "Cylinder",length=(19, "cm"),  # height = 19 cm
    radius=(30, "mm")   # radius = 30 mm (Ø = 60 mm)
    ).get_volume_and_area()
print(f"[Geometry] V = {internal_volume.item():.6g} m^3")
print(f"[Geometry] A = {surface_area.item():.6g} m^2")

# --- 2) Migrant + LDPE film ---------------------------------------------
# Example: Irganox 1076 (phenolic AO) from public DBs.
m = migrant("irganox 1076")
film = LDPE(migrant=m, l=(100, "um"))  # 100 µm
print(f"[Migrant] {m.compound} | CAS={m.CAS} | "
      f"M={m.M} g/mol | logP={m.logP}")

# --- 3) Food proxy: semi-solid, fatty sandwich ---------------------------
class Sandwich(food.realfood,food.semisolid,food.fat):name = "sandwich"
F = Sandwich(volume=internal_volume, surfacearea=surface_area,
    contacttime=(10, "days"), contacttemperature=(7, "C"),
    substance=m, simulant="ethanol")  # fatty matrix proxy

# --- 4) Equilibrium partition (Tier M2) ----------------------------------
# K_{F/P} = k_P / k_F
K_F_over_P = film.k / F.k0
print(f"[M2] K_(F/LDPE) = {K_F_over_P.item():.4g} "
       "(<1: stronger LDPE affinity)")

# --- 5) Kinetics (Tier M3): 10 d @ 7 C -----------------------------------
# SFPPy couples partition + diffusion with T,t.
kin = F.migration(film)
CF = kin.CFtarget
s = 100 * CF / m.SML
print(
 f"[M3] 10 d @ 7 C: C_F = {CF.item():.6g}, "
 f"severity = {s.item():.2f}"
 )
# Optional plots: time course + layer profiles
kin.plotCF() # kinetics
kin.plotCx() # concentration profiles

# --- 6) Step 2: +4 h @ 25 C (warming before use) -------------------------
F2 = F.copy().update( contacttime=(4, "hours"),contacttemperature=(25, "C"))
kin2 = kin >> F2 # simulation chaining with ">>"
s2 = 100 * kin2.CFtarget / m.SML # updated severity
print(f"[M3] +4 h @ 25 C: C_F = {kin2.CFtarget.item():.6g}, "
      f"severity = {s2.item():.2f}")
# Cumulative kinetics over both steps
(kin + kin2).plotCF()
