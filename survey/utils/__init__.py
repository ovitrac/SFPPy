"""
survey.utils — small shared helpers for the survey package.

Exports:
    cf_at_user_grid — interpolate a senspatankar solution's CF back
                       onto the user-requested time grid, guarding
                       against the solver's post-contact diagnostic
                       window extending sol.t beyond ttarget.

See survey/utils/cf_extract.py for the contract.
"""
from survey.utils.cf_extract import cf_at_user_grid

__all__ = ["cf_at_user_grid"]
