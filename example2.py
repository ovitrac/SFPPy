"""
===============================================================================
SFPPy Example: mass transfer from multilayer materials
===============================================================================


# Example 2: Mass Transfer Simulation in Recycled PP Bottles
------------------------------------------------------------

## Overview
This script simulates **1D mass transfer** of **toluene** from a **300 µm thick recycled PP bottle**
into a **fatty liquid food**. It also examines the **effect of adding a 30 µm PET functional barrier (FB)**.

### Simulation Steps:
1. **Define bottle geometry** (1L bottle with a body and neck).
2. **Set up polymer layers** (PP containing toluene, with/without PET FB).
3. **Define liquid food properties** (fatty liquid, storage conditions).
4. **Run mass transfer simulations**:
   - **Without a functional barrier**.
   - **With a 30 µm PET FB**.
5. **Compare migration kinetics** between both cases.
6. **Perform a systematic study** on **FB thickness** (2 µm to 60 µm).
7. **Save and print all simulation results**.

### Expected Outcomes:
✅ **Migration kinetics analysis** (with and without FB).
✅ **Optimal PET thickness selection** to reduce migration.
✅ **Ready-to-use figures and reports** for analysis.

---

@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
"""
# %% Import Dependencies
# ----------------------
import os
from patankar.loadpubchem import migrant    # Migrant online database
from patankar.geometry import Packaging3D   # 3D geometry module
import patankar.food as food                # Food contact classes
import patankar.layer as polymer            # Polymer database
from patankar.migration import senspatankar as solver # Mass transfer solver
from patankar.migration import CFSimulationContainer as store # Store results
from patankar.migration import print_figure # Printing functions
from patankar.layer import _toSI            # Convert units to SI

import matplotlib.pyplot as plt

# %% Output Folder
# ----------------
# Define the output directory to store results
outputfolder = os.path.join(os.getcwd(), "tmp")
os.makedirs(outputfolder, exist_ok=True)  # Create folder if missing

# %% Define Contact Conditions
# ----------------------------
# Store numbers with their units in a tuple
contactTemperature = (20, "degC")
contactTime = (450, "days")  # Contact duration in days

# %% Set Parameters with Uncertainty
# ----------------------------------
# Note: Concentrations have arbitrary units. Final results maintain the same units.
maxConcentration = 10  # Maximum initial concentration in the polymer (mg/kg)

# %% Define Bottle Geometry
# -------------------------
"""
Create a **1L bottle** by assembling two cylinders: **body** and **neck**.
The geometry is defined using the `Packaging3D` class.
"""
bottle = Packaging3D(
    "bottle",
    body_radius=(40, "mm"),  # 40 mm body radius
    body_height=(0.2, "m"),  # 20 cm body height
    neck_radius=(1.8, "cm"),  # 1.8 cm neck radius
    neck_height=0.05          # 5 cm neck height (default SI unit: meters)
)

# Compute internal volume (m³) and contact surface area (m²)
internalvolume, contactsurface = bottle.get_volume_and_area()

# %% Define the Migrant (Toluene)
# --------------------------------
# Retrieve chemical properties of the migrant (toluene)
surrogate = migrant("toluene")

# %% Define PP Bottle Walls
# -------------------------
# Create a **300 µm thick PP layer** containing **toluene** at the initial concentration
PPwalls_with_toluene = polymer.PP(
    l=(300, "um"),    # Thickness: 300 µm
    substance=surrogate,
    C0=maxConcentration,  # Initial concentration in PP
    T=contactTemperature  # Temperature of contact
)

# %% Define Liquid Food Properties
# --------------------------------
"""
Define a **fatty liquid food** using multiple inheritance from `food.realfood`,
`food.liquid`, and `food.fat`.
"""
class liquidFood(food.realfood, food.liquid, food.fat):
    name = "liquidFood"

# Instantiate the liquid food layer
FOODlayer = liquidFood(
    volume=internalvolume,
    surfacearea=contactsurface,
    contacttime=contactTime,
    contacttemperature=contactTemperature
)

# %% Run Mass Transfer Simulation (No Functional Barrier)
# -------------------------------------------------------
# Simulate migration for the bottle **without a functional barrier**.
ref_simulation = solver(
    PPwalls_with_toluene,  # PP walls containing toluene
    FOODlayer,             # Fatty liquid food
    name="bottle-rPP"
)

# %% Plot Results: No Functional Barrier
# --------------------------------------
# Concentration profiles at the final time
ref_simulation.plotCx()

# Migration kinetics (CF vs time) for selected time points
listoftimes = _toSI(([2,10,30,60,80,100,150,200,300,400,500,600,700], "days"))
hfig1 = ref_simulation.plotCF(t=listoftimes)

# %% Define PET Functional Barrier (FB)
# -------------------------------------
# Create a **30 µm thick PET layer** as a functional barrier
PET_functionalBarrier = polymer.gPET(
    l=(30, "um"),  # Thickness: 30 µm
    substance=surrogate,
    C0=0,  # Virgin PET (no migrant initially)
    T=contactTemperature
)

# Combine PET and PP layers (Food is on the left, PET is the outermost layer)
FBwalls_with_toluene = PET_functionalBarrier + PPwalls_with_toluene

# %% Run Simulation with Functional Barrier
# -----------------------------------------
# Simulate migration with **PET Functional Barrier**
fb_simulation = solver(
    FBwalls_with_toluene,
    FOODlayer,
    name="bottleFB-PET-rPP"
)

# %% Plot Results: With Functional Barrier
# ----------------------------------------
# Concentration profiles at final time
fb_simulation.plotCx()

# Migration kinetics (CF vs time)
hfig2 = fb_simulation.plotCF(t=listoftimes)

# %% Compare Migration: With vs Without Functional Barrier
# --------------------------------------------------------
# Store and compare both simulations
allCF = store(name="Rbottles")
allCF.add(ref_simulation, "without FB", "b")  # Reference case (blue)
allCF.add(fb_simulation, "with FB", "m")      # Functional barrier case (magenta)
hfig12 = allCF.plotCF()

# %% Systematic Study: FB Thickness (2 µm to 60 µm)
# -------------------------------------------------
"""
Evaluate migration reduction by **varying the PET barrier thickness** from **2 µm to 60 µm**
(step of 4 µm). The results are stored in `fullcomparison`.
"""
fullcomparison = store(name="fb study")
fullcomparison.add(ref_simulation, "without FB", "b")  # Add reference case

# Copy the FB configuration to modify thickness
currentfb_walls = FBwalls_with_toluene.copy()

# Iterate over thickness values from 2 µm to 60 µm in steps of 4 µm
for fb_thickness in range(2, 61, 4):
    print(f"Solving for FB = {fb_thickness} µm")

    # Update PET layer thickness (Layer 0 is the FB)
    currentfb_walls.l[0] = _toSI((fb_thickness, "um")).item()

    # Run simulation for the current FB thickness
    currentfb_simulation = solver(
        currentfb_walls,
        FOODlayer,
        name=f"bottleFB-PET-{fb_thickness}um"
    )

    # Assign a unique color for each thickness
    color = plt.cm.viridis((fb_thickness - 2) / (60 - 2))
    fullcomparison.add(currentfb_simulation, f"FB = {fb_thickness} µm", color)

# Plot migration kinetics for all FB thicknesses
hfig_all = fullcomparison.plotCF()

# %% Save and Print Figures
# -------------------------
printconfig = {"destinationfolder": outputfolder, "overwrite": True}
print_figure(hfig1, **printconfig)
print_figure(hfig2, **printconfig)
print_figure(hfig12, **printconfig)
print_figure(hfig_all, **printconfig)
