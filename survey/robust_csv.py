"""Robust CSV reader + CSV→XLSX converter — re-exported from :mod:`tools.robust_csv`.

SFPPy's core kernel (``patankar``) does not read CSV; the ``survey`` layer may.
This shim lets survey code use ``from survey.robust_csv import read_csv_robust,
to_xlsx`` while the single source of truth lives in ``tools/robust_csv.py``.

Author: Olivier Vitrac, PhD, HDR | olivier.vitrac@gmail.com | INRAE + Generative Simulation
"""
from tools.robust_csv import (  # noqa: F401
    CSVReadResult,
    read_csv_robust,
    read_csv_simple,
    to_xlsx,
)

__all__ = ["CSVReadResult", "read_csv_robust", "read_csv_simple", "to_xlsx"]


if __name__ == "__main__":  # allow `python -m survey.robust_csv <csv>...`
    import sys
    from tools.robust_csv import _main
    sys.exit(_main())
