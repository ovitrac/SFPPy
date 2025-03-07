#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
SFPPy: Approximate Correlation Between Polarity Index (P') and logP
===============================================================================

## SFPPy: Approximate Correlation Between Polarity Index (P’) and logP

### **Objective**:

Build a *qualitative* correlation between an approximate “polarity index” ($P'$) and $logP$ (octanol-water partition coefficient)
values for small organic molecules (denoted $i$) using the Flory-Huggins (FH) theory of mixtures.

$P'$ values are itended to replace the square roor of molar cohesive energies $\sqrt(E)$ in the approximation
of Flory-Huggins (FH) coefficients $i$ in a fluid or polymer phase $k$:

$$
         \chi_{i+k} \approx 0.34 \times (\sqrt{E_i} - \sqrt{E_k}) ^ 2
$$

where $E_i$ and $E_k$ are cohesive energies of $i$ and $k$, respectively. The coefficient 0.34 is empirical.

> *References*:

> - Guillaume Gillet, Olivier Vitrac, and Stéphane Desobry
>     Prediction of Solute Partition Coefficients between Polyolefins and Alcohols Using a Generalized Flory−Huggins Approach
>     *Industrial & Engineering Chemistry Research* **2009** 48 (11), 5285-5301
>     https://doi.org/10.1021/ie801141h

> - Thomas Brouwer and Boelo Schuur
>     Model Performances Evaluated for Infinite Dilution Activity Coefficients Prediction at 298.15 K
>     *Industrial & Engineering Chemistry Research* **2019** 58 (20), 8903-8914
>     https://doi.org/10.1021/acs.iecr.9b00727

Since $logP$ values include a significant entropic contribution for large substances, the effect is removed using the same FH theory.
For solute $i$, $logP = ln10 \cdot lnK_{i,o/w}$ with:

$$
    lnK_{i,o/w} \approx \gamma_{i,w} - \gamma_{i,o}
$$

At infinite dilution in $k=o,w$, we get from FH theory:

$$
        \gamma_{i,k} \approx \chi_{i+k} +1 - r_{i,k}
$$

with $r_{i,k} \approx \frac{V_i}{V_k}$ the ratio of molar volumes between $i$ and $k$

As a result:

$$
    ln10 \times logP \approx  \Delta\chi_{i}^{scale:ow} + S_i^{correction:ow}
$$

with $S_i^{correction:ow} = -\left( r_{i,w} - r_{i,o} \right)$

By replacing 0.34 by a scaling parameter $\alpha$ and by replacing cohesive energies $E$ by their polarity indices $P'$', one gets:

$$
     \Delta\chi_{i}^{scale:ow} = \chi_{i,w} - \chi_{i,o}
   \approx \alpha \times \left[ (P'_i - P'_w) ^ 2 - (P'_i - P'_o) ^ 2 \right] \\
   = \alpha \times \left[P'_i \cdot (P'_w - P'_o) + (P'_w + P'_o)(P'_w - P'_o) \\ \right] \\
   = \alpha\cdot(P'_w - P'_o)\times\left(P'_i + \beta \right) \\
   = A\cdot P'_i + B
$$

> We notice that $P'$ are related proportional to the difference of FH coefficients between
> water (w) and octanol (o) via two constants $A$ and $B$. Since $P`$ values are used as differences,
> we can arbitrarily choose $B=0$ without affecting predictions.

Commonly, $P’$ lies between $0$ (very hydrophobic substances) and $10.2$ (water = the most polar).

From these definitions, $P'$ is approximated as:

$$
    A \cdot P' \approx ln10\times logP - S_i^{correction:ow}(V_i)  \approx U\left(logP,V_i\right)
$$

> Choosing $A = 1$ leads to $\alpha = \frac{1}{(P'_w - P'_o)}$. As it will be shown and due to
> the fitting approximation, the largest amplitude leads to $\alpha\approx0.14$.

Bounding $P’$ within a finite interval and in particular the zero bound for hydrophobic substances leads several substances associated with different $logP$ (unbounded) and $V_i$ (unbounded) to have the same $P’$ values. A smooth predictor, is sought from databases by looking for a second-degree polynomial approximation of $P’$, denoted $\hat{P’}$, with $U\left(logP,V_i\right)$:

$$
	\hat{P’}(M,V) = a \cdot  \left[U\left(logP,V_i=f(M)\right)\right]^2 + b \cdot  U\left(logP,V_i=f(M)\right) + c
$$

where $V_i$ is approximated with molecular mass $M$ using the Miller’s formula:

$$
	V_i \approx 0.997 \cdot M^{1.03}
$$

With $M$ in $g \\cdot mol^{\_1} and $V_i$ in $cm^3 \cdot mol^{-1}$.

**Methodology of fitting**:

- $a$,$b$ and $c$ are fitted on a small datasets of values with $logP$ and $P’$ values reported in open databases
- The maximum amplitude is $P’$ is determined by the value for water, which is depends on the values of $logP$ reported in databases and the choices of parameters.

### **Application to the prediction of activity coefficients** 🔗 🔗 🔗 🔗 🔗 🔗 🔗

$$
\gamma_{i,k} \approx \chi_{i+k} +1 - \left(r_{i,k} - n(r_{i,k})\right)
\approx \alpha \times [\hat{P’}(logP_i,V_i) - \hat{P’}(logP_k,V_k)] ^ 2  + 1 - \left(r_{i,k} - n(r_{i,k})\right)
$$

For activity coefficients in polymers ($k=P$), $r_{i,k} - n(r_{i,k})\approx0$ since the solute is considered infinitely small compared to $P$.
In food simulants ($k=F$), $n(r_{i,k})$ represents a correction due to the presence of $n(r_{i,k}$ molecules within the spherical volume (coarse-grain approximation of FH interactions) of the solute. This number can be determined uniquely from $i+k$ pairwise radial distributions. It is larger for large flexible substances and almost $0$ for small or bulky solutes with spherical symmetry.

In the sake of a “universal” rule applicable to any substance, it is proposed to take :

$$
	n(r) = \frac{r-5}{r}
$$

> *References*:
>
> - Guillaume Gillet, Olivier Vitrac, and Stéphane Desobry
>     Prediction of Partition Coefficients of Plastic Additives between Packaging Materials and Food Simulants.
>     *Industrial & Engineering Chemistry Research* **2010** 49 (16), 7263-7280
>     https://doi.org/10.1021/ie9010595




This script
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
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from patankar.loadpubchem import migrant

# %% cache debug issues
# from patankar.loadpubchem import CompoundIndex
# dbdefault = CompoundIndex(cache_dir="cache.PubChem", index_file="pubchem_index.json")


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

# Molar volume (from its approximation with M)
V = [migrant(m).molarvolumeMiller for m in solvents]
Vwater = migrant("water").molarvolumeMiller
Voctanol = migrant("octanol").molarvolumeMiller
FHentropy_contribution = [- v * (1/Vwater - 1/Voctanol) for v in V]

# Polarity Index (P') with minor manual adjustments (denoted +)
polarity_index = [10.2,      # Water
                  5.1,     # Methanol (5.1 + 3 = 8.1)
                  4.3,       # Ethanol (4.3 + 0.7 = 5.0)
                  5.1,       # Acetone (5.1 + 0.5 = 5.6)
                  5.8,       # Acetonitrile (5.8 + 1 = 6.8)
                  3.1,       # Dichloromethane
                  2.4,   # Toluene
                  0.0        # n-Hexane
                 ]

# logP reference data
logP_values = [-1.38+0.7,  # Water
               -0.77,  # Methanol
               -0.24,  # Ethanol
               -0.21,  # Acetone
               -0.22,  # Acetonitrile
                1.25,  # Dichloromethane
                2.73,  # Toluene
                3.90   # n-Hexane
              ]

logKow_entropy = [logp * math.log(10) - c for logp, c in zip(logP_values, FHentropy_contribution)]

# % Extended Dataset for Validation
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

# Dictionary for replacements
replacements = {
    "Petroleum Ether": "Isohexane",
}

# Replace using list comprehension
ext_solvents = [replacements.get(s, s) for s in ext_solvents]


# Molar volume (from its approximation with M)
ext_V = [migrant(m).molarvolumeMiller for m in ext_solvents]
ext_FHentropy_contribution = [- v * (1/Vwater - 1/Voctanol) for v in ext_V]

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

ext_logKow = [logp * math.log(10) for logp in ext_logP_values]
ext_logKow_entropy = [logp * math.log(10) - c for logp, c in zip(ext_logP_values, ext_FHentropy_contribution)]

# Filter out solvents that appear in the tuned dataset
validation_solvents = []
validation_polarity_index = []
validation_logP_values = []
validation_FHentropy_contribution = []
validation_logKow_entropy = []

for i, solvent in enumerate(ext_solvents):
    if solvent not in solvents:
        validation_solvents.append(solvent)
        validation_polarity_index.append(ext_polarity_index[i])
        validation_logP_values.append(ext_logP_values[i])
        validation_FHentropy_contribution.append(ext_FHentropy_contribution[i])
        validation_logKow_entropy.append(ext_logKow_entropy[i])

# Create DataFrames for easy inspection
df = pd.DataFrame({
    "Solvent": solvents,
    "Polarity Index (P')": polarity_index,
    "logP": logP_values,
    "FHentropy_contribution": FHentropy_contribution,
    "logKow_entropy": logKow_entropy
}).sort_values(by="Polarity Index (P')", ascending=True)

df_validation = pd.DataFrame({
    "Solvent": validation_solvents,
    "Polarity Index (P')": validation_polarity_index,
    "logP": validation_logP_values,
    "FHentropy_contribution": validation_FHentropy_contribution,
    "logKow_entropy": validation_logKow_entropy
}).sort_values(by="Polarity Index (P')", ascending=True)


# %% Quick Visualization of the Tuned Data
plt.figure(figsize=(8, 5))
plt.scatter(polarity_index, logKow_entropy, color='blue', label="Reference Data (Tuned)")
for i, solvent in enumerate(solvents):
    plt.annotate(solvent, (polarity_index[i], logKow_entropy[i]),
                 fontsize=9, xytext=(5,5), textcoords='offset points')

plt.xlabel("Polarity Index (P')")
plt.ylabel("logKow")
plt.title("Polarity Index (P') vs. logKow (Tuned Dataset)")
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
coefficients = np.polyfit(polarity_index, logKow_entropy, 2)
quadratic_model = np.poly1d(coefficients)

# For plotting a smooth curve:
x_range = np.linspace(min(polarity_index), max(polarity_index), 100)
y_fitted = quadratic_model(x_range)
quad_eq_str = (f"logP ≈ {coefficients[0]:.4f} * (P')² "
               f"+ {coefficients[1]:.4f} * (P') "
               f"+ {coefficients[2]:.4f}")

print("Fitted coefficients (a, b, c):", coefficients)
print("Quadratic equation =>", quad_eq_str)

A,B,C = coefficients
# % Compare with Extended Validation Dataset
plt.figure(figsize=(8, 5))

# Plot the tuned set
plt.scatter(polarity_index, logKow_entropy, color='Crimson', label="Tuned (8 solvents)")
for i, solvent in enumerate(solvents):
    plt.annotate(solvent, (polarity_index[i], logKow_entropy[i]),
                 fontsize=8, xytext=(5,5), textcoords='offset points')

# Plot the validation set
plt.scatter(validation_polarity_index, validation_logKow_entropy,
            color='DeepSkyBlue', label="Validation (~35 solvents)")
for i, solvent in enumerate(validation_solvents):
    plt.annotate(solvent, (validation_polarity_index[i], validation_logKow_entropy[i]),
                 fontsize=6, xytext=(5,5), textcoords='offset points')

# Plot the fitted curve
plt.plot(x_range, y_fitted, color='Crimson', linestyle='--', label="Quadratic Fit")

plt.xlabel("Polarity Index (P')")
plt.ylabel("logKow - entropy contribution")
plt.title("Polarity Index (P') vs. logKow - entropy con\nQuadratic Model Fit & Validation")
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.show()


# %% Final implementation
import numpy as np

def polarity_index(logP=None,V=None,name=None,
                             Vw=19.588376948550433, # migrant("water").molarvolumeMiller
                             Vo=150.26143432234372, # migrant("octanol").molarvolumeMiller
                             A=0.18161296829146106, #0.14265915555430117,
                             B=-3.412678660396018, #-3.136976791114839,
                             C=14.813767205916765  # 14.52713235123395
                             ):
    """
    Computes the polarity index (P') from a given logP value and S value
    as the positive root of the quadratic fitted equation:

        E = A * (P')² + B * P' + C
        P' = (-B - sqrt(B² - 4A(C - E))) / (2A)

        with E = logP*ln(10) - S = Xw - Xo
             S = entropy contribution = - (V/Vw - V/Vo)

    Parameters:
    ----------
    logP : float, list, or np.ndarray
    V    : float, list or np.array
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
    >>> logP = migrant("anisole").logP
    >>> V = migrant("anisole").molarvolumeMiller
    >>> polarity_index(logP,V)
    8.34  # Example output

    >>> polarity_index_from_logP([-1.0, 0.5, 2.0])
    array([9.2, 4.5, 1.8])  # Example outputs
    """

    # Define valid logP range based on quadratic model limits
    Emin = C - B**2 / (4*A)  # ≈ -2.78 (theoretical minimum lnKow=E)
    Emax = C                 # 14.527 (theoretical maximum logP)
    Pmax = 10.2 # value for water

    def compute_P(logP_value,V_value):
        """Computes P' for a single logP and V value after input validation."""
        S = - (1/Vw - 1/Vo)*V_value
        E = logP_value * 2.302585092994046 - S # Chiw - Chio
        #print(f"E={E} with S={S}")
        if E < Emin:
            return Pmax  # Most polar (beyond water)
        if E > Emax:
            return 0.0  # Most hydrophobic (beyond n-Hexane)
        discriminant = B**2 - 4*A*(C - E)
        sqrt_discriminant = np.sqrt(discriminant)
        P2root = (-B - sqrt_discriminant) / (2*A)  # Always select P2
        return P2root if P2root<=Pmax else Pmax

    # for testing
    if logP is None or V is None:
        if name is None:
            raise ValueError('provide a valid name or logP, V pair values')
        from patankar.loadpubchem import migrant
        tmp = migrant(name)
        V = tmp.molarvolumeMiller
        logP = tmp.logP
    # Handle both single and multiple values efficiently
    if isinstance(logP, (list, tuple, np.ndarray)):
        return np.vectorize(compute_P)(logP,V)  # Vectorized for multiple inputs
    else:
        return compute_P(logP,V)


# % Demonstration
print("\n=== Demonstration of polarity_index_from_logP ===")
solutes = ["water","toluene","limonene", "anisole", "BHT", "Irganox 1076","ethylene","ethanol","propanol","octanol","water"]
m = migrant("BHT")
for s in solutes:
    logPref = migrant(s).logP
    P = polarity_index(name=s)
    print(f"{s}: logP={logPref}  => P'={P}") # :>5.2f


# %% Final Predictions vs Inputs
# Compute predictions
predictions = np.array([polarity_index(name=s) for s in ext_solvents])
ext_polarity_index = np.array(ext_polarity_index)  # Ensure numpy array
# Least Squares Fit without Intercept using NumPy
slope, _, _, _ = np.linalg.lstsq(ext_polarity_index.reshape(-1, 1), predictions, rcond=None)
# Generate fitted line values
x_range = np.linspace(min(ext_polarity_index), max(ext_polarity_index), 100)
y_fitted = slope[0] * x_range  # Since no intercept, y = m*x
# Plot
plt.figure(figsize=(8, 5))
plt.scatter(ext_polarity_index, predictions, color='Crimson', label="Tuned (8 solvents)")
# Annotate solvents
for i, solvent in enumerate(ext_solvents):
    plt.annotate(solvent, (ext_polarity_index[i], predictions[i]),
                 fontsize=10, xytext=(5, 5), textcoords='offset points')
# 45-degree reference line
plt.plot(x_range, x_range, linestyle="--", color="black", label="45° line (y=x)")
# Fitted regression line (passing through origin)
plt.plot(x_range, y_fitted, linestyle="-", color="blue", label="Fitted line (no intercept)")
# Labels and styling
plt.xlabel("Polarity Index (P') - reference")
plt.ylabel("Polarity Index (P') - predicted")
plt.title("Predicted vs Reference P' Values")
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.show()

# %%
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
