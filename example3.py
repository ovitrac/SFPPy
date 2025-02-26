#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
SFPPy Example: mass transfer from chained steps and complex scenarios
===============================================================================


Example 3: Advanced Migration Simulation with Variants
------------------------------------------------------

Scenario: Mass Transfer in a Multilayer Packaging System
--------------------------------------------------------
This example simulates the migration of **limonene** from a recycled **polypropylene (PP)** core
to food, using a **trilayer packaging (ABA structure)** where **A = PET** and **B = PP**.
The study follows three consecutive storage and processing steps:

1. Storage at ambient temperature (20°C) for 4 months
2. Hot-filling at 90°C with a fatty food
3. Storage at 30°C for 6 months

Limonene is initially present in **B at 200 mg/kg**. The layer thicknesses are:
- PET (A): **20 µm**
- PP (B): **500 µm**

Variants Studied:
- **Variant 1**: Replace limonene with toluene
- **Variant 2**: Reduce B’s thickness from 20 µm to 10 µm for both limonene and toluene
- **Variant 3**: Combination of Variant 1 and Variant 2

All four conditions (Reference, Variant 1, Variant 2, and Variant 3) are compared in a single
figure, with results exported to **Excel and CSV** for further analysis.

--------------------------------------------------------
Exploring SFPPy’s Database
--------------------------------------------------------
This example allows interaction with **SFPPy’s extensive databases**, including:
- **3D geometries** (predefined packaging shapes)
- **Chemical substances** (millions via PubChem)
- **Food contact conditions** (various storage and processing scenarios)
- **Polymers and materials** (tens of predefined types)

--------------------------------------------------------
Key Features and Pythonic Syntax
--------------------------------------------------------

1. **Compact and Concise Code Structure**
Example 3 demonstrates **chained simulations** with **minimal coding**, contrasting with
Examples 1 and 2. The approach **stores results within objects**, enabling streamlined operations.

2. **Operator `>>` for Seamless Processing**
The `>>` operator automates **property propagation** and **mass transfer simulations**:
- `mypackaging >> myfood` → Updates **food geometry**
- `food1 >> food2` → Synchronizes food properties
- `myfood >> mymaterial` → Transfers **temperature**
- `mymaterial >> myfood` → Simulates **mass transfer**

Example: `mypackaging >> myfood >> mymaterial >> myfood` performs:
1. **Update food geometry** to match packaging
2. **Transfer temperature** to the material
3. **Simulate mass transfer** from material to food

For multiple contact conditions:
`mypackaging >> foodcontact2 >> foodcontact1 >> mymaterial >> foodcontact1 >> foodcontact2`
- Adjusts food geometry (`foodcontact1`, `foodcontact2`)
- Transfers **temperature** (`foodcontact1 >> mymaterial`)
- Simulates mass transfer and propagates results
- Final results stored in **foodcontact1.lastsimulation, foodcontact2.lastsimulation**

To explicitly save results: `result = ... step n >> step n+1`
Dynamic updates:
`foodcontact1 >> mymaterial.update(substance="new migrant", l=new_thickness) >> foodcontact1`

3. **Operator `+` for Combining Results**
- **Material assembly:** `ABA = A + B + A` (trilayer structure)
- **Merging simulation results:** `fullsolution = foodcontact1.lastsimulation + foodcontact2.lastsimulation`

4. **Comprehensive Result Output**
- **Visualization:** Results can be plotted directly
- **Exporting:** Save to **Excel** (`.xlsx`) and **CSV** (`.csv`)

--------------------------------------------------------
Summary
--------------------------------------------------------
Example 3 demonstrates **efficient chaining of migration simulations** using intuitive operators
(`>>` for sequential processing, `+` for merging results). SFPPy’s powerful database and
Pythonic syntax enable complex simulations with minimal effort.


@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
"""

# %% Output folder
# -----------------
# Define the output directory to store results.
import os
outputfolder = os.path.join(os.getcwd(), "tmp")  # Full path
os.makedirs(outputfolder, exist_ok=True)  # Create folder if missing

# %% Build the geometry
# ---------------------
"""
Identify and select a suitable 3D packaging geometry.
The function `help_geometry()` lists available shapes.

Example: `box_container` is a rectangular prism with one open face.
Its volume and surface area are computed as:
    Volume = length × width × height
    Surface Area = 2(lw + lh + wh) - lw  (removing the open top)
"""
# Import the geometry module
from patankar.geometry import Packaging3D, help_geometry

# Display available geometries and required parameters
help_geometry()

# Define the container as an open rectangular prism
container = Packaging3D('box_container', height=(8, "cm"), width=(10, "cm"), length=(19, "cm"))

# %% Define the migrant
# ---------------------
# Retrieve the properties of the migrating substance.
from patankar.loadpubchem import migrant

# Search for limonene in PubChem and retrieve chemical properties.
m = migrant("limonene")

# %% Identify the materials
# -------------------------
"""
List available polymer/material options for packaging.
Use `help_layer()` to display predefined materials.

| Class Name | Type     | Material                 | Code |
|------------|---------|---------------------------|------|
| PP         | polymer | Isotactic Polypropylene   | PP   |
| gPET       | polymer | Glassy PET                | PET  |
"""
# Import material properties
from patankar.layer import help_layer, layer

# Display available polymer options
help_layer()

# Import required polymer classes (Polypropylene, PET)
from patankar.layer import gPET, PP

# %% Build the materials
# ----------------------
"""
Define and configure the layers for the multilayer packaging.
Each layer has properties such as:
    - Thickness (`l`)
    - Initial concentration (`C0`)
    - Migrant type (`substance`, `migrant`, etc.)

Property synonyms:
| Parameter  | Synonyms                  |
|------------|---------------------------|
| substance  | molecule, solute, migrant |
| C0         | Cp0, CP0                  |
| l          | lp, lP                    |
"""
# Define the PET (A) and PP (B) layers
A = gPET(l=(20, "um"), migrant=m, C0=0)       # 20 µm PET with no initial migrant
B = PP(l=(0.5, "mm"), migrant=m, CP0=200)    # 500 µm PP with 200 mg/kg limonene

# Assemble the multilayer structure: ABA configuration
ABA = A + B + A

# Display basic information (default temperature: 40°C)
print("\nOur ABA technology\n", repr(ABA))

# %% Identify storage and contact conditions
# -----------------------------------------
"""
Define different storage and processing conditions.
Use `help_food()` to display available predefined food-contact scenarios.

| Class Name  | Description              | Level  | Inheritance    | Parameters        |
|-------------|--------------------------|--------|----------------|-------------------|
| fat        | Fat contact               | root   |chemicalaffinity| k0                |
| liquid     | Liquid food               | root   | texture        | h                 |
| ambient    | Ambient storage           | contact| realcontact    | contacttemperature|
| hotfilled  | Hot-filling               | contact| realcontact    | contacttemperature|
| stacked    | Stacked storage           | user   | setoff         | h                 |
"""
from patankar.food import help_food
help_food()

# Import relevant food-contact conditions
from patankar.food import realfood, liquid, fat, ambient, hotfilled, stacked

# Define three food-contact conditions
class contact1(stacked, ambient):
    name = "1:setoff"
    contacttemperature = (20, "degC")
    contacttime = (3, "months")

class contact2(hotfilled, realfood, liquid, fat):
    name = "2:hotfilling"

class contact3(ambient, realfood, liquid, fat):
    name = "3:storage"
    contacttime = (6, "months")

# Instantiate food-contact conditions
medium1 = contact1(contacttime=(4, "months"))  # Adjust time
medium2 = contact2()
medium3 = contact3()

# %% Wrap the food with its container-packaging
# --------------------------------------------
"""
Use the `>>` operator to propagate the packaging properties (volume, surface area)
to all food-contact conditions.

Example:
    container >> medium1 >> medium2 >> medium3
"""
container >> medium1 >> medium2 >> medium3

# %% Simulate Step 1
# ------------------
"""
Chain operations using `>>`:
1. Transfer temperature from `medium1` to `ABA`
2. Simulate mass transfer from `ABA` to `medium1`
"""
medium1 >> ABA >> medium1

# Verify simulation results
medium1.lastsimulation.plotCx()

# %% Simulate all steps (Step 1 → Step 2 → Step 3)
# -----------------------------------------------
"""
Continue simulation by chaining:
    medium1 >> ABA >> medium1 >> medium2 >> medium3

All results are stored inside:
    medium1.lastsimulation, medium2.lastsimulation, medium3.lastsimulation
"""
medium1 >> ABA >> medium1 >> medium2 >> medium3

# Plot migration kinetics
medium1.lastsimulation.plotCx()
medium2.lastsimulation.plotCF()
medium3.lastsimulation.plotCF()

# %% Combine all profiles into a single kinetic curve
# ---------------------------------------------------
"""
Use `+` operator to merge kinetic profiles.
"""
sol123 = medium1.lastsimulation + medium2.lastsimulation + medium3.lastsimulation
sol123.plotCF()

# %% Variant 1: Replace limonene with toluene
# -------------------------------------------
"""
Repeat all steps using toluene instead of limonene.
"""
m2 = migrant("toluene")  # Retrieve new migrant

# Restart the simulation pipeline with updated migrant
medium1 >> ABA.update(solute=m2) >> medium1 >> medium2 >> medium3

# Store results
sol123_variant1 = medium1.lastsimulation + medium2.lastsimulation + medium3.lastsimulation
sol123_variant1.plotCF()

# %% Variant 2: Reduce thickness of outer PET layers
# --------------------------------------------------
"""
Halve the thickness of the first and last layers (A) while keeping PP constant.
"""
refthickness = ABA.l.copy()
newthickness = refthickness
newthickness[[0, -1]] /= 2  # Reduce thickness of A layers (index 0=first layer in contact with food, index -1= last layer)

# Restart simulation with modified structure
medium1 >> ABA.copy(l=newthickness, migrant=m) >> medium1 >> medium2 >> medium3
sol123_variant2 = medium1.lastsimulation + medium2.lastsimulation + medium3.lastsimulation

# %% Variant 3: Combine Variant 1 and Variant 2
# ---------------------------------------------
"""
Use both modifications:
- Toluene as the migrant
- Reduced PET layer thickness
"""
medium1 >> ABA.copy(l=newthickness, migrant=m2) >> medium1 >> medium2 >> medium3
sol123_variant3 = medium1.lastsimulation + medium2.lastsimulation + medium3.lastsimulation

# %% Compare Reference and Variants
# ---------------------------------
"""
Store and compare all solutions using `CFSimulationContainer`.
"""
from patankar.migration import CFSimulationContainer as store

# Initialize collection
collection = store(name="container with FB")

# Add results to collection
collection.add(sol123, "Limonene-FB = 20 µm (ref)", "dodgerblue", linewidth=3)
collection.add(sol123_variant1, "Toluene-FB = 20 µm (v1)", "orangered", linewidth=3)
collection.add(sol123_variant2, "Limonene-FB = 10 µm (v2)", "deepskyblue", linewidth=2)
collection.add(sol123_variant3, "Toluene-FB = 10 µm (v3)", "tomato", linewidth=2)

# Plot comparative migration kinetics
collection_fig = collection.plotCF()

# %% Export Results (PNG, PDF, Excel, CSV)
# ----------------------------------------
"""
Save simulation results for further analysis.
"""
printconfig = {"destinationfolder": outputfolder, "overwrite": True}
collection_fig.print(**printconfig)

# Export data
collection.save_as_excel(filename="example3.xlsx", destinationfolder=outputfolder, overwrite=True)
collection.save_as_csv(filename="example3.csv", destinationfolder=outputfolder, overwrite=True)
