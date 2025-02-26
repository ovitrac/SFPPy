"""
	PATANKAR package
	Initialization of the Patankar package (prevent multiple instances of generic libraries)

	@author: olivier.vitrac@agroparistech.fr

	History
	2021-01-23 first generation

"""

# Package Dependencies (put in __init__.py)
# =========================================
import numpy as np
from copy import deepcopy as duplicate
from patankar.pint import UnitRegistry as SIbase
from patankar.struct import struct
SI = SIbase()
qSI = SI.Quantity


# all variable
# ============
__all__ = ["struct", "layer","senspatankar"]
