"""
Regression tests for the SFPPylite (Pyodide) memory guard of patankar.migration.

Under Pyodide (`_LITE_` True), a SensPatankarResult keeps no restart
interpolator: the PCHIP used by resumeat() stores four coefficient arrays the
size of the whole Cx(t) history and exhausts WebAssembly memory. The tests check
that the guard is inert on desktop, that the lite branch drops the interpolator
and disables resumeat() with an explicit message, and that the simulated values,
resume() and the restart file are unchanged.

@project: SFPPy - Safe Food Packaging in Python
@author: Olivier Vitrac
@license: MIT
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.interpolate import PchipInterpolator

import patankar.migration as migration
from patankar.food import ethanol
from patankar.layer import layer
from patankar.migration import senspatankar


def _solve():
    """Bilayer, fresh solution over 25 days (as in the migration.py demo)."""
    A = layer(layername="layer A", k=2, C0=0, D=1e-16)
    B = layer(layername="layer B")
    return senspatankar(A + B, ethanol(), t=(25, "days"))


@pytest.fixture
def lite(monkeypatch):
    """Run the code path taken under Pyodide."""
    monkeypatch.setattr(migration, "_LITE_", True)


def test_guard_inert_on_desktop():
    assert migration._LITE_ is False
    sol = _solve()
    assert isinstance(sol._restart_Cxi_interp, PchipInterpolator)
    sol2 = sol.resumeat((10, "days"))
    assert np.all(np.isfinite(sol2.CF))


def test_lite_keeps_no_interpolator(lite):
    sol = _solve()
    assert sol._restart_Cxi_interp is None
    with pytest.raises(ValueError, match="SFPPylite"):
        sol.resumeat((10, "days"))


def test_lite_results_and_resume_unchanged(monkeypatch):
    ref = _solve()
    ref2 = ref.resume((40, "days"))
    monkeypatch.setattr(migration, "_LITE_", True)
    sol = _solve()
    assert_allclose(sol.CF, ref.CF, rtol=0, atol=0)
    assert_allclose(sol.Cx, ref.Cx, rtol=0, atol=0)
    sol2 = sol.resume((40, "days"))
    assert_allclose(sol2.CF, ref2.CF, rtol=1e-12)
