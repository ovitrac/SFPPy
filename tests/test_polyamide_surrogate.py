"""
Regression test for the monomer surrogate of the polyamides (patankar.layer PA6, PA66).

PA6 (-NH-(CH2)5-CO-) and PA66 (-NH-(CH2)6-NH-CO-(CH2)4-CO-) carry the same amide density,
one amide per six carbons, so both use n-hexanamide. Their Henry-like coefficients then differ
only through crystallinity: k is scaled by 1/(1 - c). Before SFPPy 1.9.3, PA66 used adipamide
(two amides per six carbons), which inflated k for apolar solutes about 50-fold.

@project: SFPPy - Safe Food Packaging in Python
@author: Olivier Vitrac
@license: MIT
"""

import numpy as np
from numpy.testing import assert_allclose

import patankar.layer as polymer
from patankar.loadpubchem import migrant


def test_pa6_pa66_same_surrogate():
    assert polymer.PA6().chemicalsubstance == polymer.PA66().chemicalsubstance == "n-hexanamide"


def test_pa66_differs_from_pa6_by_crystallinity_only():
    m = migrant("limonene")
    pa6 = polymer.PA6(substance=m, T=(40, "degC"))
    pa66 = polymer.PA66(substance=m, T=(40, "degC"))
    c6, c66 = pa6.crystallinity_history[0], pa66.crystallinity_history[0]
    ratio = float(np.asarray(pa66.k).flat[0]) / float(np.asarray(pa6.k).flat[0])
    assert_allclose(ratio, (1 - c6) / (1 - c66), rtol=1e-10)
