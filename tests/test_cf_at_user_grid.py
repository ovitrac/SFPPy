"""
survey.utils.cf_at_user_grid helper invariants.

Runs a minimal senspatankar call, verifies:

  (a) sol.CF[-1] is NOT at the user's ttarget (the post-contact window
      makes it diverge).
  (b) cf_at_user_grid(sol, [ttarget])[0] agrees with sol.CFtarget.
  (c) cf_at_user_grid(sol, user_grid) is monotone-non-decreasing on
      intermediate points when the underlying physics is monotone.

These are fast integration tests (~1-2 s) — the solver call itself is
tiny (simple monolayer, short t).

@project: SFPPy — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def _minimal_monolayer_sol(t_grid, contacttime=None):
    """
    Run a quick senspatankar on a tiny LDPE monolayer + ethanol.

    When `contacttime` is None, it is taken as max(t_grid) — the
    "no extension" case. When contacttime > max(t_grid), the solver
    extends sol.t with a post-contact diagnostic window and sol.t[-1]
    sits at (or past) contacttime. This is what triggers the
    sol.CF[-1] ≠ sol.CFtarget divergence.
    """
    from patankar.layer import layer
    from patankar.food import foodphysics
    from patankar.migration import senspatankar
    t_grid = np.asarray(t_grid, dtype=float)
    if contacttime is None:
        contacttime = float(t_grid.max())
    lay = layer(l=1e-4, D=1e-13, k=1.0, C0=1.0, T=20.0)
    food = foodphysics(
        k=1.0, k0=1.0, h=1e-7,
        surfacearea=1.0, volume=10.0,
        contacttime=float(contacttime),
        contacttemperature=20.0, CF0=0.0,
    )
    return senspatankar(lay, food, t=t_grid, autotime=False)


def test_cf_at_user_grid_matches_CFtarget_when_contacttime_matches_t():
    """
    When `food.contacttime == max(t)`, `sol.ttarget == max(t)` and
    `cf_at_user_grid(sol, [max(t)])` agrees with `sol.CFtarget`.
    """
    from survey.utils import cf_at_user_grid
    t_user = np.array([0.0, 1e-15])
    sol = _minimal_monolayer_sol(t_user)  # contacttime defaults to max(t)

    cf_helper = float(cf_at_user_grid(sol, [t_user[-1]])[0])
    cf_target = float(sol.CFtarget)

    scale = max(abs(cf_target), 1e-15)
    assert abs(cf_helper - cf_target) <= 1e-10 * scale, \
        f"helper {cf_helper} vs CFtarget {cf_target}"


def test_cf_at_user_grid_documents_the_extension_pitfall():
    """
    DOCUMENTS the pitfall: when `food.contacttime` > max(user `t`), the
    solver extends `sol.t` with a post-contact diagnostic window.
    `sol.CF[-1]` then sits at or past `contacttime`, NOT at max(t).
    `sol.CFtarget` interpolates at `sol.ttarget = contacttime` — which
    also is NOT max(t) under this call pattern.

    The safe reads under `contacttime > max(t)`:
        np.interp(t_user_last, sol.t, sol.CF)      # explicit interp
        cf_at_user_grid(sol, [t_user_last])[0]     # canonical helper

    Callers that want `CFtarget == CF(max(t))` must pass
    `contacttime = max(t)` to foodphysics (this is what the survey workers do
    in its per-iteration food construction).
    """
    t_user = np.array([0.0, 1e-15])
    # Force contacttime >> max(t) — triggers the extension.
    sol = _minimal_monolayer_sol(t_user, contacttime=1.5)

    cf_last = float(np.asarray(sol.CF)[-1])
    cf_target = float(sol.CFtarget)
    cf_at_user_last = float(
        np.interp(t_user[-1], np.asarray(sol.t), np.asarray(sol.CF)))

    # ttarget = contacttime = 1.5, NOT max(t) = 1e-15
    assert float(np.asarray(sol.ttarget)[0]) == 1.5

    # CFtarget reads CF at contacttime (1.5), very different from CF at
    # user's t=1e-15. This is the trap.
    assert cf_target > 1000.0 * cf_at_user_last, (
        f"CFtarget = {cf_target} should be >>CF at t=1e-15 "
        f"= {cf_at_user_last} when contacttime >> max(t). "
        f"If this assertion fires, senspatankar's ttarget semantics "
        f"may have changed — review survey/utils/cf_extract.py.")

    # sol.CF[-1] sits at the extended endpoint (1.8 in our case),
    # even further from the user's requested t=1e-15 than CFtarget.
    assert cf_last >= cf_target, (
        f"sol.CF[-1] = {cf_last} should be >= CFtarget = {cf_target} "
        f"under the extension; if not, revisit the solver's "
        f"post-contact diagnostic window.")


def test_cf_at_user_grid_preserves_monotonicity():
    """
    CF(t) for a monotone-releasing monolayer must be non-decreasing on
    the user's grid. Helper preserves that invariant.
    """
    from survey.utils import cf_at_user_grid
    t_user = np.linspace(0.0, 5e4, 20)  # seconds
    sol = _minimal_monolayer_sol(t_user)
    cf = cf_at_user_grid(sol, t_user)
    dcf = np.diff(cf)
    assert np.all(dcf >= -1e-15), \
        f"CF not monotone-non-decreasing on user grid: dCF min = {dcf.min()}"


def test_cf_at_user_grid_shape_mismatch_raises():
    """Helper validates input shapes."""
    from survey.utils import cf_at_user_grid
    class _FakeSol:
        t = np.array([0.0, 1.0])
        CF = np.array([0.0, 1.0, 2.0])  # mismatched
        CFtarget = 1.0
    with pytest.raises(ValueError, match="shape mismatch"):
        cf_at_user_grid(_FakeSol(), [0.5])


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
