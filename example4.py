#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
SFPPy Example: fitting an arbirary combination of parameters D[i] and k[i]
from experiments using layerLinks and dynamic simulations.
===============================================================================

Example 4: Sensitivity Analysis and Kinetic Data Fitting
---------------------------------------------------------

# **Example4.py: Estimation of D and k from Concentration Kinetics in Food Simulants**
## **Overview**
This script illustrates how **diffusivity (D)** and **partition coefficients (k)** can be estimated from **experimental concentration kinetics** in **food simulants**. The method relies on **layerLink objects**, which allow adjusting simulation parameters dynamically.

The workflow involves:
    1. **Define a monolayer material (P) and a food simulant (F).**
    2. **Enforce a reference partition coefficient (kF = 1) for F**, allowing kP to be interpreted as the **partition coefficient between F and P**.
    3. **Perform a migration simulation** to establish a baseline.
    4. **Generate pseudo-experimental data** with controlled noise levels.
    5. **Evaluate the squared error (`d2 = R - E`)** between simulation (`R`) and experiment (`E`).
    6. **Conduct a sensitivity analysis** by varying D and k to assess their impact.
    7. **Optimize D and k simultaneously** to recover their original values.
    8. **Analyze the goodness of fit and risk of overfitting.**

## **Key Concepts**
- **`layerLink` for parameter linking:** Enables modifying **D** and **k** externally without modifying the base layer object.
- **Migration simulation (`F.migration(P)`)**: Computes mass transfer from P into F.
- **Pseudo-experimental data (`R.pseudoexperiment()`)**: Creates synthetic datasets for validation.
- **Fitting (`R.fit(E)`)**: Automatically adjusts **D** and **k** to best match experimental data.

## **Expected Outcomes**
- **Validation of layerLink objects** for parameter optimization.
- **Comparison of simulated vs. experimental data** to evaluate fit quality.
- **Determination of optimal D and k values** for given experimental conditions.

--------------------------------------------------------
Summary
--------------------------------------------------------
This script demonstrates **an advanced approach for parameter estimation** in food safety analysis.

@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
"""

# %% Import Dependencies
from patankar.layer import layer, layerLink
from patankar.food import foodlayer

# %% Define a Monolayer Material (P)
P = layer(
    l=(100, "um"),
    D=(1e-10, "cm**2/s"),  # Diffusivity in cm²/s (internally converted to SI)
    C0=1000,  # Initial concentration (arbitrary units)
    k=0.1  # Partition coefficient
)

# %% Define a Generic Contact Medium (F)
F = foodlayer(
    contacttime=(10, "days"),
    volume=(1, "L"),
    surfacearea=(6, "dm**2"),
    h=(1e-6, "m/s"),  # Mass transfer coefficient
    CF0=0,  # Initial concentration in the food
    k=1  # Fixed kF, allowing kP to represent kF/P
)

# %% Create and Attach `layerLink` Objects for D and k
Dreference, kreference = P.D, P.k  # Store original values for comparison
D = layerLink("D", indices=0, values=Dreference)  # Create link for D
k = layerLink("k", indices=0, values=kreference)  # Create link for k
P.Dlink = D  # Attach links to P
P.klink = k

# %% Run Initial Migration Simulation
R = F.migration(P)  # Perform mass transfer simulation
R.plotCF()  # Plot the reference concentration profile in food

# %% Generate Pseudo-Experimental Data
E = R.pseudoexperiment(npoints=30, std_relative=0.01)  # Simulated experimental dataset
E.plotCF()

# %% Compute Initial Squared Distance (d2)
d2 = R - E  # Compute squared error function
d2_original = d2()  # Evaluate initial deviation from experiment
print(f"Initial squared distance (E-R)**2: {d2_original}")

# %% Sensitivity Analysis: Explore the Impact of D and k Variations
R.comparison.update(0, label="**Reference**", color="black")
R.comparison.plotCF()

# Set iteration count and color mapping
niterations = 10
cmap = R.comparison.jet(niterations)

# Vary D and k systematically
for i in range(1, niterations + 1):
    D[0] = D[0] / 1.1  # Decrease D by 10%
    k[0] = k[0] * 1.1  # Increase k by 10%

    # Format LaTeX labels for D and k
    name = f"{P.Dlatex()[0]}, {P.klatex()[0]}"

    # Rerun simulation with updated parameters
    R.rerun(name=name, color=cmap[i-1])

    # Compute new squared distance
    d2new = d2()
    print(f"[{i}/{niterations}]: Distance variation = {100 * d2new / d2_original - 100:.2f}%")

# Add pseudo-experimental data to comparison
R.comparison.add(E, label="Pseudo Experiment", discrete=True)
R.comparison.plotCF()  # Final visualization

# %% Optimize D and k to Recover the Original Values
d2beforeOptim = d2()  # Distance before fitting
resfit = R.fit(E)  # Perform parameter fitting
d2afterOptim = d2()  # Distance after fitting

# Print optimization results
print(f"BEFORE OPTIMIZATION: Distance (E-R)**2 = {d2beforeOptim}")
print(f"AFTER OPTIMIZATION: Distance (E-R)**2 = {d2afterOptim}")
overfitting_flag = "[Overfitting detected]" if d2afterOptim < d2_original else ""
print(f"Variation with original = {100 * d2afterOptim / d2_original - 100:.4g}% {overfitting_flag}\n")
print(f"Fitted D = {D.values} [m²/s] vs. Original D = {Dreference} [m²/s]")
print(f"Fitted k = {k.values} [a.u.] vs. Original k = {kreference} [a.u.]")

# Plot optimized results
R.comparison.plotCF()
