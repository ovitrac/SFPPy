"""
Pytest configuration and shared fixtures for SFPPy tests.

@project: SFPPy - Safe Food Packaging in Python
@author: Olivier Vitrac
@license: MIT
"""

import sys
from pathlib import Path

import pytest
import numpy as np

# Ensure SFPPy root is in path
_sfppy_root = Path(__file__).resolve().parent.parent
if str(_sfppy_root) not in sys.path:
    sys.path.insert(0, str(_sfppy_root))


# =============================================================================
# Fixtures for layer tests
# =============================================================================

@pytest.fixture
def ldpe_layer():
    """Create a basic LDPE layer for testing."""
    from patankar.layer import LDPE
    return LDPE(l=(100, "um"), C0=1000, T=(25, "degC"))


@pytest.fixture
def hdpe_layer():
    """Create a basic HDPE layer for testing."""
    from patankar.layer import HDPE
    return HDPE(l=(50, "um"), C0=500, T=(40, "degC"))


@pytest.fixture
def gpet_layer():
    """Create a basic gPET layer for testing."""
    from patankar.layer import gPET
    return gPET(l=(12, "um"), C0=0, T=(25, "degC"))


# =============================================================================
# Fixtures for food tests
# =============================================================================

@pytest.fixture
def ethanol_food():
    """Create ethanol food simulant."""
    from patankar.food import ethanol
    return ethanol(
        contacttime=(10, "days"),
        contacttemperature=(25, "degC")
    )


@pytest.fixture
def water_food():
    """Create water food simulant."""
    from patankar.food import water
    return water(
        contacttime=(10, "days"),
        contacttemperature=(25, "degC")
    )


# =============================================================================
# Fixtures for migrant tests
# =============================================================================

@pytest.fixture
def limonene_migrant():
    """Create limonene migrant."""
    from patankar.loadpubchem import migrant
    return migrant("limonene")


@pytest.fixture
def toluene_migrant():
    """Create toluene migrant."""
    from patankar.loadpubchem import migrant
    return migrant("toluene")


# =============================================================================
# Numerical tolerance fixtures
# =============================================================================

@pytest.fixture
def rtol():
    """Relative tolerance for floating point comparisons."""
    return 1e-6


@pytest.fixture
def atol():
    """Absolute tolerance for floating point comparisons."""
    return 1e-12
