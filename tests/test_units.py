"""
Unit conversion tests for SFPPy.

Tests the _toSI() and check_units() functions which are critical
for correct physical calculations throughout the codebase.

@project: SFPPy - Safe Food Packaging in Python
@author: Olivier Vitrac
@license: MIT
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from patankar.layer import _toSI, check_units


class TestToSI:
    """Tests for the _toSI() unit conversion function."""

    # =========================================================================
    # Length conversions
    # =========================================================================

    def test_micrometers_to_meters(self):
        """100 µm = 1e-4 m"""
        result = _toSI((100, "um"))
        assert_allclose(result, 1e-4, rtol=1e-10)

    def test_millimeters_to_meters(self):
        """1 mm = 1e-3 m"""
        result = _toSI((1, "mm"))
        assert_allclose(result, 1e-3, rtol=1e-10)

    def test_centimeters_to_meters(self):
        """10 cm = 0.1 m"""
        result = _toSI((10, "cm"))
        assert_allclose(result, 0.1, rtol=1e-10)

    def test_nanometers_to_meters(self):
        """500 nm = 5e-7 m"""
        result = _toSI((500, "nm"))
        assert_allclose(result, 5e-7, rtol=1e-10)

    # =========================================================================
    # Time conversions
    # =========================================================================

    def test_days_to_seconds(self):
        """10 days = 864000 s"""
        result = _toSI((10, "days"))
        assert_allclose(result, 864000, rtol=1e-10)

    def test_hours_to_seconds(self):
        """24 hours = 86400 s"""
        result = _toSI((24, "hours"))
        assert_allclose(result, 86400, rtol=1e-10)

    def test_minutes_to_seconds(self):
        """60 minutes = 3600 s"""
        result = _toSI((60, "minutes"))
        assert_allclose(result, 3600, rtol=1e-10)

    def test_months_to_seconds(self):
        """1 month ~ 30.44 days = 2629746 s (average)"""
        result = _toSI((1, "months"))
        # Month is approximately 30.44 days
        assert result > 2.5e6 and result < 2.7e6

    def test_years_to_seconds(self):
        """1 year = 31557600 s (365.25 days)"""
        result = _toSI((1, "years"))
        assert_allclose(result, 31557600, rtol=1e-3)

    # =========================================================================
    # Diffusivity conversions (area/time)
    # =========================================================================

    def test_cm2_per_s_to_m2_per_s(self):
        """1 cm²/s = 1e-4 m²/s"""
        result = _toSI((1, "cm**2/s"))
        assert_allclose(result, 1e-4, rtol=1e-10)

    def test_um2_per_s_to_m2_per_s(self):
        """1 µm²/s = 1e-12 m²/s"""
        result = _toSI((1, "um**2/s"))
        assert_allclose(result, 1e-12, rtol=1e-10)

    # =========================================================================
    # Temperature conversions
    # =========================================================================

    def test_celsius_passthrough(self):
        """25°C should remain 25 (default output is degC)"""
        result = _toSI((25, "degC"))
        assert_allclose(result, 25, rtol=1e-10)

    def test_kelvin_to_celsius(self):
        """298.15 K = 25°C"""
        result = _toSI((298.15, "K"))
        # Note: check_units returns degC by default for temperatures
        assert_allclose(result, 25, rtol=1e-3)

    def test_freezing_point_kelvin(self):
        """273.15 K = 0°C"""
        result = _toSI((273.15, "K"))
        assert_allclose(result, 0, atol=0.01)

    # =========================================================================
    # Volume conversions
    # =========================================================================

    def test_liters_to_m3(self):
        """1 L = 1e-3 m³"""
        result = _toSI((1, "L"))
        assert_allclose(result, 1e-3, rtol=1e-10)

    def test_milliliters_to_m3(self):
        """1000 mL = 1e-3 m³"""
        result = _toSI((1000, "mL"))
        assert_allclose(result, 1e-3, rtol=1e-10)

    # =========================================================================
    # Area conversions
    # =========================================================================

    def test_cm2_to_m2(self):
        """100 cm² = 0.01 m²"""
        result = _toSI((100, "cm**2"))
        assert_allclose(result, 0.01, rtol=1e-10)

    # =========================================================================
    # Array inputs
    # =========================================================================

    def test_list_input(self):
        """Lists should be converted element-wise."""
        result = _toSI(([1, 2, 3], "mm"))
        expected = np.array([1e-3, 2e-3, 3e-3])
        assert_allclose(result.flatten(), expected, rtol=1e-10)

    def test_numpy_array_input(self):
        """NumPy arrays should be converted element-wise."""
        result = _toSI((np.array([100, 200]), "um"))
        expected = np.array([1e-4, 2e-4])
        assert_allclose(result.flatten(), expected, rtol=1e-10)

    # =========================================================================
    # Error handling
    # =========================================================================

    def test_invalid_tuple_raises(self):
        """Non-tuple input should raise ValueError."""
        with pytest.raises(ValueError):
            _toSI(100)

    def test_invalid_tuple_length_raises(self):
        """Tuple with wrong length should raise ValueError."""
        with pytest.raises(ValueError):
            _toSI((100, "um", "extra"))

    def test_invalid_unit_type_raises(self):
        """Non-string unit should raise ValueError."""
        with pytest.raises(ValueError):
            _toSI((100, 123))


class TestCheckUnits:
    """Tests for the check_units() function."""

    def test_tuple_input(self):
        """Tuple input should be parsed correctly."""
        value, units = check_units((100, "um"))
        assert_allclose(value, 1e-4, rtol=1e-10)

    def test_numpy_array_passthrough(self):
        """NumPy arrays should pass through unchanged."""
        arr = np.array([1, 2, 3])
        result, units = check_units(arr)
        assert np.array_equal(result, arr)

    def test_none_passthrough(self):
        """None should pass through unchanged."""
        result, units = check_units(None)
        assert result is None

    def test_temperature_celsius_to_kelvin(self):
        """Temperature conversion from degC to K."""
        value, units = check_units(
            (25, "degC"),
            ExpectedUnits="K"
        )
        assert_allclose(value, 298.15, rtol=1e-3)
        assert units == "K"

    def test_temperature_kelvin_to_celsius(self):
        """Temperature conversion from K to degC."""
        value, units = check_units(
            (298.15, "K"),
            ExpectedUnits="degC"
        )
        assert_allclose(value, 25, rtol=1e-3)
        assert units == "degC"


class TestUnitConsistency:
    """Tests to verify unit consistency across typical use cases."""

    def test_layer_thickness_typical_range(self):
        """Typical packaging thicknesses should convert correctly."""
        # Thin film: 12 µm
        thin = _toSI((12, "um"))
        assert 1e-5 < thin < 1e-4

        # Thick bottle: 500 µm
        thick = _toSI((500, "um"))
        assert 4e-4 < thick < 6e-4

        # Very thick: 2 mm
        very_thick = _toSI((2, "mm"))
        assert_allclose(very_thick, 2e-3, rtol=1e-10)

    def test_contact_time_typical_range(self):
        """Typical contact times should convert correctly."""
        # Short contact: 2 hours
        short = _toSI((2, "hours"))
        assert_allclose(short, 7200, rtol=1e-10)

        # Standard test: 10 days
        standard = _toSI((10, "days"))
        assert_allclose(standard, 864000, rtol=1e-10)

        # Long storage: 6 months
        long = _toSI((6, "months"))
        assert long > 1.5e7  # > 6 months in seconds

    def test_diffusivity_typical_range(self):
        """Typical diffusivity values should convert correctly."""
        # Fast diffusion in rubber: 1e-8 cm²/s
        fast = _toSI((1e-8, "cm**2/s"))
        assert_allclose(fast, 1e-12, rtol=1e-10)

        # Slow diffusion in PET: 1e-14 cm²/s
        slow = _toSI((1e-14, "cm**2/s"))
        assert_allclose(slow, 1e-18, rtol=1e-10)
