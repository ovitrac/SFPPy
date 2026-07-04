"""
survey/utils/cf_extract.py — safe extraction of CF from a
senspatankar result onto a user-specified time grid.

Background
----------
senspatankar appends a post-contact diagnostic window to sol.t (used
internally to compute PRT and PR_effective). As a result:

    sol.t[-1] > ttarget       for every non-degenerate call
    sol.CF[-1]                is CF at sol.t[-1], NOT at ttarget

Consumers who read sol.CF[-1] expecting "CF at the user's requested
final time" will silently get a value near the equilibrium CF
(because the diagnostic window sits well past contacttime). This was
the root cause of a two-step extraction bug fixed during survey-engine
validation.

Rule
----
    cf_at_ttarget           = sol.CFtarget          (canonical, interp_CF)
    cf_at_arbitrary_grid(g) = cf_at_user_grid(sol, g)    (this module)

Never use sol.CF[-1] or sol.CF[i] to read "CF at the user's
time[-1] / time[i]" — the solver's extension breaks that assumption.

@project: SFPPy — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""
from __future__ import annotations

from typing import Sequence, Union

import numpy as np


def cf_at_user_grid(
    sol,
    user_grid: Union[Sequence[float], np.ndarray],
) -> np.ndarray:
    """
    Interpolate a senspatankar solution's CF onto a user-supplied grid.

    Parameters
    ----------
    sol : SensPatankarResult
        Result of a senspatankar() call. Must expose `.t` and `.CF`
        (and `.CFtarget` as a sanity anchor).
    user_grid : sequence of float
        Time points at which CF is requested. Values outside
        [sol.t.min(), sol.t.max()] are clamped to the boundary values
        (same as np.interp default).

    Returns
    -------
    np.ndarray
        CF values at each user_grid point, shape (len(user_grid),).

    Notes
    -----
    - Safe regardless of whether the solver extended sol.t beyond
      the user's ttarget.
    - For a single-point read at ttarget, prefer sol.CFtarget
      directly (no interpolation overhead).
    - This function does NOT verify that user_grid is a subset of
      sol.t — np.interp silently interpolates or clamps. Callers that
      need strict grid alignment should compare user_grid against
      sol.t themselves.
    """
    t = np.asarray(sol.t, dtype=float).reshape(-1)
    cf = np.asarray(sol.CF, dtype=float).reshape(-1)
    grid = np.asarray(user_grid, dtype=float).reshape(-1)
    if t.shape != cf.shape:
        raise ValueError(
            f"sol.t and sol.CF shape mismatch: {t.shape} vs {cf.shape}"
        )
    # np.interp requires monotonically non-decreasing xp.
    order = np.argsort(t)
    return np.interp(grid, t[order], cf[order])
