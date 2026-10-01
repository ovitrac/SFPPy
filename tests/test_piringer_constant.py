"""
Regression test for the Piringer diffusivity model (patankar.property.Dpiringer).

The reference equation (JRC practical guidelines, EUR 27529 EN, 2015; Zhu, Welle &
Vitrac, Soft Matter 15(42):8912-8932, 2019, eq. 43) reads, in m2/s,

    D = exp(A''_P - tau/T - 0.1351 M^(2/3) + 0.003 M - 10454/T)

SFPPy used 0.135 before 1.9.

@project: SFPPy - Safe Food Packaging in Python
@author: Olivier Vitrac
@license: MIT
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from patankar.property import Dpiringer


def reference(App, tau, M, T):
    TK = T + 273.15
    return np.exp(App - tau / TK - 0.1351 * M ** (2.0 / 3.0) + 0.003 * M - 10454.0 / TK)


@pytest.mark.parametrize("polymer, App, tau", [("LDPE", 11.5, 0), ("HDPE", 14.5, 1577)])
@pytest.mark.parametrize("M, T", [(100.0, 40.0), (500.0, 60.0), (1200.0, 20.0)])
def test_piringer_reference_equation(polymer, App, tau, M, T):
    assert_allclose(Dpiringer.evaluate(polymer=polymer, M=M, T=T), reference(App, tau, M, T), rtol=1e-12)
    assert_allclose(Dpiringer(polymer=polymer).eval(M=M, T=T), reference(App, tau, M, T), rtol=1e-12)
