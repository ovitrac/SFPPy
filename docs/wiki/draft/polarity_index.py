#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
SFPPy: Approximate Correlation Between Polarity Index (P') and logP
===============================================================================
**Objective**:
    Build a *qualitative* correlation between an approximate "polarity index"
    (P') and logP (octanol-water partition coefficient) values. This script
    demonstrates how to:

    1) Fit a *quadratic* model using a small, "tuned" dataset that covers a
       range of polarities from n-Hexane to Water.
    2) Validate it (roughly) with an extended dataset of ~35 solvents.
    3) Provide a function to invert the model and estimate P' from logP values,
       for classification or for ranking in subsequent Henry-like or partition
       coefficient models in SFPPy.

**Disclaimer**:
    - This correlation is *qualitative* and based on an intentionally *small*
      dataset, hand-tuned for a few representative solvents (non-polar to very
      polar).
    - Do **not** expect high accuracy. The goal is to produce a rough scale (P')
      that loosely tracks how "polar" or "non-polar" a compound might be.
    - The script demonstrates the approach rather than guaranteeing a universal
      model.

**References & Data**:
    - Data for Polarity Index (P') adapted from various solvent polarity tables:
      LSU Macromolecular Resources, HPLC Solvent Guides, and additional web
      resources.
    - logP values drawn from typical reference tables (PubChem, various
      compiled datasheets).

**Motivation**:
    We wish to provide a "fast" way to guess an ordering of polarity from logP,
    which is widely available for many chemicals (e.g., in PubChem). P' is then
    used within SFPPy to set or guess Henry-like coefficients or Flory–Huggins
    parameters in patankar.layer and patankar.food modules.

----------------------------------------------------------------------------
Created on Fri Feb 28 12:01:56 2025
@author: olivi (community)
----------------------------------------------------------------------------
"""

# %% Import libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# %% Tuned Reference Dataset
"""
This small set of 8 solvents is carefully chosen to span from n-Hexane (non-polar)
to Water (very polar). Polarity Index (P') values are lightly 'tweaked' to get a
smooth progression. The logP values come from known literature data.
"""

solvents = [
    "Water", "Methanol", "Ethanol", "Acetone", "Acetonitrile",
    "Dichloromethane", "Toluene", "n-Hexane"
]

# Polarity Index (P') with minor manual adjustments (denoted +)
polarity_index = [10.2,          # Water
                  5.1+3,         # Methanol (5.1 + 3 = 8.1)
                  4.3+0.7,       # Ethanol (4.3 + 0.7 = 5.0)
                  5.1+0.5,       # Acetone (5.1 + 0.5 = 5.6)
                  5.8+1,         # Acetonitrile (5.8 + 1 = 6.8)
                  3.1,           # Dichloromethane
                  2.4,           # Toluene
                  0.0            # n-Hexane
                 ]

# logP reference data
logP_values = [-1.38,  # Water
               -0.77,  # Methanol
               -0.24,  # Ethanol
               -0.21,  # Acetone
               -0.22,  # Acetonitrile
                1.25,  # Dichloromethane
                2.73,  # Toluene
                3.90   # n-Hexane
              ]

# %% Extended Dataset for Validation
"""
This larger set (~35 solvents) is used to see if the small dataset's fit
extends (qualitatively) to other solvents. We exclude the ones already in
the small set to avoid double-counting.
"""

ext_solvents = [
    "Pentane", "1,1,2-Trichlorotrifluoroethane", "Cyclopentane", "Heptane", "Hexane",
    "Iso-Octane", "Petroleum Ether", "Cyclohexane", "n-Butyl Chloride", "Toluene",
    "Methyl t-Butyl Ether", "o-Xylene", "Chlorobenzene", "o-Dichlorobenzene", "Ethyl Ether",
    "Dichloromethane", "Ethylene Dichloride", "n-Butyl Alcohol", "Isopropyl Alcohol",
    "n-Butyl Acetate", "Isobutyl Alcohol", "Methyl Isoamyl Ketone", "n-Propyl Alcohol",
    "Tetrahydrofuran", "Chloroform", "Methyl Isobutyl Ketone", "Ethyl Acetate",
    "Methyl n-Propyl Ketone", "Methyl Ethyl Ketone", "1,4-Dioxane", "Acetone", "Methanol",
    "Pyridine", "2-Methoxyethanol", "Acetonitrile", "Propylene Carbonate", "N,N-Dimethylformamide",
    "Dimethyl Acetamide", "N-Methylpyrrolidone", "Dimethyl Sulfoxide", "Water"
]

ext_polarity_index = [
    0.0, 0.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 1.0, 2.4, 2.5, 2.5, 2.7, 2.7, 2.8,
    3.1, 3.5, 3.9, 3.9, 4.0, 4.0, 4.0, 4.0, 4.0, 4.1, 4.2, 4.4, 4.5, 4.7, 4.8,
    5.1, 5.1, 5.3, 5.5, 5.8, 6.1, 6.4, 6.5, 6.7, 7.2, 10.2
]

ext_logP_values = [
    3.39, 4.30, 3.20, 4.66, 3.90, 4.50, 3.50, 3.44, 2.70, 2.73, 1.20, 3.12, 2.84, 3.38, 0.83,
    1.25, 1.48, 0.88, 0.05, 1.82, 0.79, 1.98, 0.25, 0.46, 1.97, 1.31, 0.73, 1.50, 0.29, -0.27,
    -0.24, -0.77, 0.65, -0.77, -0.22, -0.41, -1.01, -0.77, -0.38, -1.35, -1.38
]

# Filter out solvents that appear in the tuned dataset
validation_solvents = []
validation_polarity_index = []
validation_logP_values = []

for i, solvent in enumerate(ext_solvents):
    if solvent not in solvents:
        validation_solvents.append(solvent)
        validation_polarity_index.append(ext_polarity_index[i])
        validation_logP_values.append(ext_logP_values[i])

# Create DataFrames for easy inspection
df = pd.DataFrame({
    "Solvent": solvents,
    "Polarity Index (P')": polarity_index,
    "logP": logP_values
}).sort_values(by="Polarity Index (P')", ascending=True)

df_validation = pd.DataFrame({
    "Solvent": validation_solvents,
    "Polarity Index (P')": validation_polarity_index,
    "logP": validation_logP_values
}).sort_values(by="Polarity Index (P')", ascending=True)


# %% Quick Visualization of the Tuned Data
plt.figure(figsize=(8, 5))
plt.scatter(polarity_index, logP_values, color='blue', label="Reference Data (Tuned)")
for i, solvent in enumerate(solvents):
    plt.annotate(solvent, (polarity_index[i], logP_values[i]),
                 fontsize=9, xytext=(5,5), textcoords='offset points')

plt.xlabel("Polarity Index (P')")
plt.ylabel("logP")
plt.title("Polarity Index (P') vs. logP (Tuned Dataset)")
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)  # Zero line for logP
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.show()


# %% Fit a Quadratic Model
"""
We assume a model:

    logP = a * (P')^2 + b * P' + c

We'll use np.polyfit() with degree=2. Then we compare it with
the extended validation set.
"""
coefficients = np.polyfit(polarity_index, logP_values, 2)
quadratic_model = np.poly1d(coefficients)

# For plotting a smooth curve:
x_range = np.linspace(min(polarity_index), max(polarity_index), 100)
y_fitted = quadratic_model(x_range)
quad_eq_str = (f"logP ≈ {coefficients[0]:.4f} * (P')² "
               f"+ {coefficients[1]:.4f} * (P') "
               f"+ {coefficients[2]:.4f}")

print("Fitted coefficients (a, b, c):", coefficients)
print("Quadratic equation =>", quad_eq_str)


# %% Compare with Extended Validation Dataset
plt.figure(figsize=(8, 5))

# Plot the tuned set
plt.scatter(polarity_index, logP_values, color='Crimson', label="Tuned (8 solvents)")
for i, solvent in enumerate(solvents):
    plt.annotate(solvent, (polarity_index[i], logP_values[i]),
                 fontsize=8, xytext=(5,5), textcoords='offset points')

# Plot the validation set
plt.scatter(validation_polarity_index, validation_logP_values,
            color='DeepSkyBlue', label="Validation (~35 solvents)")
for i, solvent in enumerate(validation_solvents):
    plt.annotate(solvent, (validation_polarity_index[i], validation_logP_values[i]),
                 fontsize=6, xytext=(5,5), textcoords='offset points')

# Plot the fitted curve
plt.plot(x_range, y_fitted, color='Crimson', linestyle='--', label="Quadratic Fit")

plt.xlabel("Polarity Index (P')")
plt.ylabel("logP")
plt.title("Polarity Index (P') vs. logP\nQuadratic Model Fit & Validation")
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.show()


# %% Final implementation
import numpy as np

def polarity_index_from_logP(logP,
                             A=0.04426556625879231,
                             B=-0.9796466111259537,
                             C=4.086222276379432):
    """
    Computes the polarity index (P') from a given logP value
    as the positive root of the quadratic fitted equation:

        logP = A * (P')² + B * P' + C
        P' = (-B - sqrt(B² - 4A(C - logP))) / (2A)

    Parameters:
    ----------
    logP : float, list, or np.ndarray
        The logP value(s) for which to compute the polarity index P'.

    Returns:
    -------
    float or np.ndarray
        The calculated polarity index P' corresponding to the input logP.
        If logP is out of the valid range, returns:
        - 10.2 for very polar solvents (beyond water)
        - 0 for extremely hydrophobic solvents (beyond n-Hexane)

    Example Usage:
    -------------
    >>> polarity_index_from_logP(-0.5)
    8.34  # Example output

    >>> polarity_index_from_logP([-1.0, 0.5, 2.0])
    array([9.2, 4.5, 1.8])  # Example outputs
    """

    # Define valid logP range based on quadratic model limits
    logPmin = C - B**2 / (4*A)  # ≈ -1.334 (theoretical minimum logP)
    logPmax = C                  # 4.086 (theoretical maximum logP)
    Pmax = 10.2 # value for water

    def compute_P(logP_value):
        """Computes P' for a single logP value after input validation."""
        if logP_value < logPmin:
            return Pmax  # Most polar (beyond water)
        if logP_value > logPmax:
            return 0.0  # Most hydrophobic (beyond n-Hexane)

        discriminant = B**2 - 4*A*(C - logP_value)
        sqrt_discriminant = np.sqrt(discriminant)
        P2root = (-B - sqrt_discriminant) / (2*A)  # Always select P2
        return P2root if P2root<=Pmax else Pmax

    # Handle both single and multiple values efficiently
    if isinstance(logP, (list, tuple, np.ndarray)):
        return np.vectorize(compute_P)(logP)  # Vectorized for multiple inputs
    else:
        return compute_P(logP)


# %% Demonstration
print("\n=== Demonstration of polarity_index_from_logP ===")
demo_values = [-1.5, -0.5, 0.5, 5.0]
results = polarity_index_from_logP(demo_values)
for lv, r in zip(demo_values, results):
    print(f"logP={lv:>5.2f} => P'={r:>5.2f}")


# %% Final Plot - Extended
plt.figure(figsize=(8, 5))
plt.scatter(ext_polarity_index, ext_logP_values, color='Teal', label="Extended Solvents")
for i, solvent in enumerate(ext_solvents):
    plt.annotate(solvent, (ext_polarity_index[i], ext_logP_values[i]),
                 fontsize=7, xytext=(5,5), textcoords='offset points')

logP_range = np.linspace(-4, 6, 1000)
p_estimated = polarity_index_from_logP(logP_range)
plt.plot(p_estimated, logP_range, color='Crimson', linestyle='-.',
         label="Inverse Quadratic Model (P'(logP))")

plt.xlabel("Polarity Index (P')")
plt.ylabel("logP")
plt.title("Approximate Mapping: P' <-> logP (Extended View)")
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.show()

"""
===============================================================================
Summary & Notes:
----------------
1) We fitted a simple quadratic (logP ~ a(P')² + bP' + c) to a small set
   of 8 solvents. The shape is roughly correct for a wide range of polarities.

2) The extended dataset (~35 solvents) is used for a rough *visual* validation.
   There's no guarantee of accuracy, especially for more exotic or borderline
   solvents.

3) The final function `polarity_index_from_logP()` is the main deliverable.
   - It clamps P' between 0.0 (extremely non-polar) and 10.2 (extremely polar).
   - Perfect for a quick rank-based approach or a first-guess at how "polar"
     an unknown solvent might be relative to water, methanol, or hexane.

4) Additional disclaimers:
   - The model is purely *empirical*.
   - Use at your own risk for approximate classification or quick bounding.

===============================================================================
"""
