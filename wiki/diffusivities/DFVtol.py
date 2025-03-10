import numpy as np
import matplotlib.pyplot as plt

# This script shows the prediction of diffusivities of toluene using hole Free-Volume theory
# at different temperatures in different polymers (LDPE, PMMA, PS, PVAc, PET).
# gPET corresponds to dry PET (unplasticized) and wPET corresponds to plasticized PET.
#
# Data have been extracted from the script figure8_revSep2019_2025.m used to generate
# Figure 9 for the reference mentioned below.
#
# REFERENCE
# Zhu Y., Welle, F. and Vitrac O. A blob model to parameterize polymer hole free volumes and solute diffusion",
# *Soft Matter* **2019**, 15(42), 8912-8932. DOI: https://doi.org/10.1039/C9SM01556F
#
# ABSTRACT
# Solute diffusion in solid polymers has tremendous applications in packaging,
# reservoir, and biomedical technologies but remains poorly understood. Diffusion
# of non-entangled linear solutes with chemically identical patterns (blobs) deviates
# dramatically in polymers in the solid-state (αlin > 1, Macromolecules 2013, 46, 874)
# from their behaviors in the molten state (αlin = 1, Macromolecules, 2007, 40, 3970).
# This work uses the scale invariance of the diffusivities, D, of linear probes
# D(N·M_blob + M_anchor,T,Tg) = N^(-αlin(T,Tg)) * D(M_blob + M_anchor,T,Tg) comprising
# N identical blobs of mass M_blob and possibly one different terminal pattern (anchor of
# mass M_anchor) to evaluate the amounts of hole-free volume in seven polymers (aliphatic,
# semi-aromatic and aromatic) over a broad range of temperatures (−70 K ≤ T − Tg ≤ 160 K).
# The new parameterization of the concept of hole-free volumes opens the application of
# the free-volume theory (FVT) developed by Vrentas and Duda to practically any polymer,
# regardless of the availability of free-volume parameters. The quality of the estimations
# was tested with various probes including n-alkanes, 1-alcohols, n-alkyl acetates, and
# n-alkylbenzene. The effects of enthalpic and entropic effects of the blobs and the anchor
# were analyzed and quantified. Blind validation of the reformulated FVT was tested
# successfully by predicting from first principles the diffusivities of water and toluene
# in amorphous polyethylene terephthalate from 4 °C to 180 °C and in various other polymers.
# The new blob model would open the rational design of additives with controlled diffusivities
# in thermoplastics.

# Constants
R = 8.31
T0K = 273.15  # K

# Polymer data (Tg in K) stored in a dictionary.
data = {
    'LDPE': {'Tg': 148.15, 'D0': 1.87e-08, 'xi': 0.615,  'ref': 3, 'Ka': 144, 'Kb': 40, 'E': 0, 'r': 0.5},
    'PMMA': {'Tg': 381.15, 'D0': 1.87e-08, 'xi': 0.56,   'ref': 2, 'Ka': 252, 'Kb': 65, 'E': 0, 'r': 0.5},
    'PS':   {'Tg': 373.15, 'D0': 4.8e-08,  'xi': 0.584,  'ref': 2, 'Ka': 144, 'Kb': 40, 'E': 0, 'r': 0.5},
    'PVAc': {'Tg': 305.15, 'D0': 1.87e-08, 'xi': 0.86,   'ref': 4, 'Ka': 142, 'Kb': 40, 'E': 0, 'r': 0.5},
    'gPET': {'Tg': 349.15, 'D0': 1.0205e-08, 'xi': 0.6761, 'ref': 5, 'Ka': 252, 'Kb': 65, 'E': 0, 'r': 0.6153},
    'wPET': {'Tg': 316.15, 'D0': 1.02046e-08, 'xi': 0.6761, 'ref': 5, 'Ka': 252, 'Kb': 65, 'E': 0, 'r': 0.277734},
}

# References list (for documentation)
references = [
    'Vrentas and Vrentas, 1994',
    'Zielinski and Duda, 1992',
    'Lutzow et al., 1999',
    'Hong, 1995',
    # toluene in PET
    'Welle,2008',
    'Pennarun et al., 2004',
    'Welle,2013',
    'our work (permeation)',
    'our work (sorption)',
]

# Main function that computes the predicted diffusivity DscalingPlike
def DscalingPlike(T, polymer):
    """
    Compute the predicted diffusivity using the scaling law.

    Parameters:
        T : float
            Temperature in Kelvin.
        polymer : str
            Polymer name (must be one of the keys in the data dictionary).

    Returns:
        float: Predicted diffusivity.
    """
    # Helper function to lookup a property value from the data dictionary.
    def lookup(P, prop):
        return data[P][prop]

    # Define alpha for T >= Tg.
    def alpha(T, Tg, Ka, Kb):
        return 1 + Ka / (T - Tg + Kb)

    # Define alphag for T < Tg.
    def alphag(T, Tg, Ka, Kb, r):
        return 1 + Ka / (r * (T - Tg) + Kb)

    deltaT = 2  # (K) sharpness of the transition at Tg

    # Heaviside-like function using tanh.
    def H(T_val, Tg_val):
        return 0.5 * (1 + np.tanh(4 / deltaT * (T_val - Tg_val)))

    # Composite alpha function that smoothly transitions between alpha and alphag.
    def alphaT(T, Tg, Ka, Kb, r):
        # Note: H(lookup(polymer, 'Tg'), T) computes H(Tg,T) and H(T, lookup(polymer, 'Tg')) computes H(T,Tg)
        return H(lookup(polymer, 'Tg'), T) * alphag(T, Tg, Ka, Kb, r) + H(T, lookup(polymer, 'Tg')) * alpha(T, Tg, Ka, Kb)

    betalin = 1
    # Function to compute Plike (see publication, Plike holds polymer effects)
    def Plike(T, polymer):
        Tg = lookup(polymer, 'Tg')
        Ka = lookup(polymer, 'Ka')
        Kb = lookup(polymer, 'Kb')
        r = lookup(polymer, 'r')
        return (alphaT(T, Tg, Ka, Kb, r) + betalin) / 0.24

    # Compute the predicted diffusivity.
    D0 = lookup(polymer, 'D0')
    E = lookup(polymer, 'E')
    xi = lookup(polymer, 'xi')

    return D0 * np.exp(-E / (R * T)) * np.exp(-xi * Plike(T, polymer))


# Test data (similar to a DataFrame) stored as a list of dictionaries.
test_data = [
    {'polymer': 'PS',    'T': 451.15, 'Dexp': 4e-11},
    {'polymer': 'PS',    'T': 433.15, 'Dexp': 1.5e-11},
    {'polymer': 'PS',    'T': 413.15, 'Dexp': 2e-12},
    {'polymer': 'PS',    'T': 383.15, 'Dexp': 1.5e-13},
    {'polymer': 'PMMA',  'T': 433.15, 'Dexp': 9e-13},
    {'polymer': 'PMMA',  'T': 413.15, 'Dexp': 2e-13},
    {'polymer': 'PMMA',  'T': 403.15, 'Dexp': 8e-14},
    {'polymer': 'LDPE',  'T': 343.15, 'Dexp': 7e-11},
    {'polymer': 'PVAc',  'T': 383.15, 'Dexp': 2e-12},
    {'polymer': 'PVAc',  'T': 353.15, 'Dexp': 1e-13},
    {'polymer': 'PVAc',  'T': 313.15, 'Dexp': 8e-16},
    {'polymer': 'PS',    'T': 388.15, 'Dexp': 7.2e-13},
    {'polymer': 'PS',    'T': 451.15, 'Dexp': 4e-11},
    {'polymer': 'PS',    'T': 433.15, 'Dexp': 1.5e-11},
    {'polymer': 'PS',    'T': 413.15, 'Dexp': 2e-12},
    {'polymer': 'PS',    'T': 383.15, 'Dexp': 1.5e-13},
    {'polymer': 'PS',    'T': 388.15, 'Dexp': 7.2e-13},
    {'polymer': 'PMMA',  'T': 433.15, 'Dexp': 9e-13},
    {'polymer': 'PMMA',  'T': 413.15, 'Dexp': 2e-13},
    {'polymer': 'PMMA',  'T': 403.15, 'Dexp': 8e-14},
    {'polymer': 'wPET',  'T': 313.15, 'Dexp': 4.1e-16},
    {'polymer': 'gPET',  'T': 313.15, 'Dexp': 2e-18},
    {'polymer': 'gPET',  'T': 313.15, 'Dexp': 3.8e-18},
    {'polymer': 'gPET',  'T': 333.15, 'Dexp': 2.6e-17},
    {'polymer': 'wPET',  'T': 313.15, 'Dexp': 6.3e-17},
    {'polymer': 'wPET',  'T': 313.15, 'Dexp': 2.1e-16},
    {'polymer': 'wPET',  'T': 313.15, 'Dexp': 4.8e-16},
    {'polymer': 'wPET',  'T': 323.15, 'Dexp': 4.1e-16},
    {'polymer': 'wPET',  'T': 323.15, 'Dexp': 6.1e-16},
    {'polymer': 'gPET',  'T': 373.15, 'Dexp': 3.26e-15},
    {'polymer': 'gPET',  'T': 453.15, 'Dexp': 2.78e-12},
    {'polymer': 'gPET',  'T': 434.15, 'Dexp': 1.17e-12},
    {'polymer': 'gPET',  'T': 414.15, 'Dexp': 1.9e-13},
    {'polymer': 'gPET',  'T': 394.15, 'Dexp': 2.47e-14},
]

# Compute predicted diffusivity for each test data point.
for entry in test_data:
    entry['Dpred'] = DscalingPlike(entry['T'], entry['polymer'])

# Plot experimental vs predicted D values.
plt.figure(figsize=(8, 6))

# Define a color map for polymers.
colors = {
    'LDPE': 'blue',
    'PMMA': 'green',
    'PS': 'red',
    'PVAc': 'orange',
    'gPET': 'purple',
    'wPET': 'brown'
}

# Plot points for each polymer with its own color.
for polymer in colors.keys():
    x = [entry['Dexp'] for entry in test_data if entry['polymer'] == polymer]
    y = [entry['Dpred'] for entry in test_data if entry['polymer'] == polymer]
    if x:  # Only plot if there are data points.
        plt.scatter(x, y, label=polymer, color=colors[polymer], edgecolors='k')

# Plot the reference line y = x.
x_vals = np.logspace(-20, -10, 100)
plt.plot(x_vals, x_vals, 'k-', label='y = x')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('experimental D values (m$^2$/s)')
plt.ylabel('predicted D values (m$^2$/s)')
plt.title('DFV vs. experimental D values')
plt.legend()
plt.show()
