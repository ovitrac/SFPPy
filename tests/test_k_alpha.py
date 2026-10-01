"""
Regression test for the scaling constant alpha of the Henry-like model (patankar.property.kFHP).

chi = alpha * (P'_i - P'_k)^2, with alpha = 0.162331 recalibrated on eight reference solvents.
Before SFPPy 1.9.1, the migrant template passed the ancient value 0.14, which overrode the
model default: layer.k was computed with alpha = 0.14.

@project: SFPPy - Safe Food Packaging in Python
@author: Olivier Vitrac
@license: MIT
"""

import inspect

import numpy as np
from numpy.testing import assert_allclose

import patankar.layer as polymer
from patankar.loadpubchem import migrant
from patankar.property import gFHP, kFHP

ALPHA = 0.162331


def test_alpha_defaults():
    for model in (kFHP, gFHP):
        assert inspect.signature(model.evaluate).parameters["alpha"].default == ALPHA
    assert migrant("limonene").ktemplate["alpha"] == ALPHA


def test_layer_k_uses_recalibrated_alpha():
    m = migrant("limonene")
    lay = polymer.PP(substance=m, T=(40, "degC"))
    mono = migrant(lay.chemicalsubstance)
    tpl = dict(m.ktemplate, Pk=mono.polarityindex, Vk=mono.molarvolumeMiller,
               crystallinity=lay.crystallinity_history[0], porosity=lay.porosity_history[0])
    k_layer = float(np.asarray(lay.k).flat[0])
    assert_allclose(k_layer, float(np.asarray(kFHP.evaluate(**dict(tpl, alpha=ALPHA))).flat[0]), rtol=1e-12)
    assert not np.isclose(k_layer, float(np.asarray(kFHP.evaluate(**dict(tpl, alpha=0.14))).flat[0]))
