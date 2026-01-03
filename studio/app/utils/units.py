"""
Unit Conversion Utilities for SFPPy Studio

Handles conversion between common units used in migration modeling.

@author: SFPPy Studio
@license: MIT
"""

from typing import Dict, Tuple, Optional
import math

# =============================================================================
# Unit Conversion Factors to SI
# =============================================================================

# Length units -> meters
LENGTH_TO_SI: Dict[str, float] = {
    "m": 1.0,
    "meter": 1.0,
    "meters": 1.0,
    "cm": 1e-2,
    "mm": 1e-3,
    "um": 1e-6,
    "µm": 1e-6,
    "micron": 1e-6,
    "microns": 1e-6,
    "nm": 1e-9,
    "in": 0.0254,
    "inch": 0.0254,
    "inches": 0.0254,
}

# Time units -> seconds
TIME_TO_SI: Dict[str, float] = {
    "s": 1.0,
    "sec": 1.0,
    "second": 1.0,
    "seconds": 1.0,
    "min": 60.0,
    "minute": 60.0,
    "minutes": 60.0,
    "h": 3600.0,
    "hr": 3600.0,
    "hour": 3600.0,
    "hours": 3600.0,
    "d": 86400.0,
    "day": 86400.0,
    "days": 86400.0,
    "week": 604800.0,
    "weeks": 604800.0,
    "month": 2592000.0,  # 30 days
    "months": 2592000.0,
    "year": 31536000.0,  # 365 days
    "years": 31536000.0,
}

# Temperature units (special handling required)
TEMPERATURE_UNITS = {"K", "degC", "C", "degF", "F", "celsius", "fahrenheit", "kelvin"}

# Mass units -> kg
MASS_TO_SI: Dict[str, float] = {
    "kg": 1.0,
    "g": 1e-3,
    "mg": 1e-6,
    "ug": 1e-9,
    "µg": 1e-9,
    "lb": 0.453592,
}

# Concentration units (dimensionless ratios)
CONCENTRATION_UNITS = ["mg/kg", "ppm", "g/kg", "µg/kg", "ppb", "%", "mg/L", "g/L"]


# =============================================================================
# Conversion Functions
# =============================================================================

def convert_length(value: float, from_unit: str, to_unit: str = "m") -> float:
    """
    Convert length between units.

    Args:
        value: Numeric value to convert
        from_unit: Source unit (e.g., "um", "mm", "cm")
        to_unit: Target unit (default "m" for SI)

    Returns:
        Converted value

    Raises:
        ValueError: If units are not recognized
    """
    from_unit = from_unit.lower().strip()
    to_unit = to_unit.lower().strip()

    if from_unit not in LENGTH_TO_SI:
        raise ValueError(f"Unknown length unit: {from_unit}")
    if to_unit not in LENGTH_TO_SI:
        raise ValueError(f"Unknown length unit: {to_unit}")

    # Convert to SI (meters), then to target unit
    si_value = value * LENGTH_TO_SI[from_unit]
    return si_value / LENGTH_TO_SI[to_unit]


def convert_time(value: float, from_unit: str, to_unit: str = "s") -> float:
    """
    Convert time between units.

    Args:
        value: Numeric value to convert
        from_unit: Source unit (e.g., "days", "hours", "months")
        to_unit: Target unit (default "s" for SI)

    Returns:
        Converted value

    Raises:
        ValueError: If units are not recognized
    """
    from_unit = from_unit.lower().strip()
    to_unit = to_unit.lower().strip()

    if from_unit not in TIME_TO_SI:
        raise ValueError(f"Unknown time unit: {from_unit}")
    if to_unit not in TIME_TO_SI:
        raise ValueError(f"Unknown time unit: {to_unit}")

    # Convert to SI (seconds), then to target unit
    si_value = value * TIME_TO_SI[from_unit]
    return si_value / TIME_TO_SI[to_unit]


def convert_temperature(value: float, from_unit: str, to_unit: str = "K") -> float:
    """
    Convert temperature between units.

    Args:
        value: Numeric value to convert
        from_unit: Source unit ("K", "degC", "C", "degF", "F")
        to_unit: Target unit (default "K" for SI)

    Returns:
        Converted value

    Raises:
        ValueError: If units are not recognized
    """
    from_unit = from_unit.strip()
    to_unit = to_unit.strip()

    # Normalize unit names
    from_norm = from_unit.lower().replace("deg", "").replace("celsius", "c").replace("fahrenheit", "f").replace("kelvin", "k")
    to_norm = to_unit.lower().replace("deg", "").replace("celsius", "c").replace("fahrenheit", "f").replace("kelvin", "k")

    # Convert to Kelvin first
    if from_norm in ("c", "degc"):
        kelvin = value + 273.15
    elif from_norm in ("f", "degf"):
        kelvin = (value - 32) * 5/9 + 273.15
    elif from_norm == "k":
        kelvin = value
    else:
        raise ValueError(f"Unknown temperature unit: {from_unit}")

    # Convert from Kelvin to target
    if to_norm in ("c", "degc"):
        return kelvin - 273.15
    elif to_norm in ("f", "degf"):
        return (kelvin - 273.15) * 9/5 + 32
    elif to_norm == "k":
        return kelvin
    else:
        raise ValueError(f"Unknown temperature unit: {to_unit}")


def convert_to_si(value: float, unit: str) -> float:
    """
    Convert any value to its SI equivalent.

    Automatically detects unit type (length, time, temperature).

    Args:
        value: Numeric value
        unit: Unit string

    Returns:
        Value in SI units (m, s, K)
    """
    unit_lower = unit.lower().strip()

    # Check length units
    if unit_lower in LENGTH_TO_SI:
        return value * LENGTH_TO_SI[unit_lower]

    # Check time units
    if unit_lower in TIME_TO_SI:
        return value * TIME_TO_SI[unit_lower]

    # Check temperature units
    if unit.strip() in TEMPERATURE_UNITS or unit_lower.replace("deg", "") in ("c", "f", "k"):
        return convert_temperature(value, unit, "K")

    raise ValueError(f"Unknown unit: {unit}")


def format_value_with_unit(value: float, unit: str, precision: int = 4) -> str:
    """
    Format a value with its unit for display.

    Uses scientific notation for very small or large values.

    Args:
        value: Numeric value
        unit: Unit string
        precision: Number of significant figures

    Returns:
        Formatted string like "1.234e-14 m²/s"
    """
    if abs(value) < 1e-3 or abs(value) > 1e4:
        return f"{value:.{precision}e} {unit}"
    else:
        return f"{value:.{precision}g} {unit}"


def get_display_unit(si_unit: str, quantity: str) -> Tuple[str, float]:
    """
    Get a human-friendly display unit and conversion factor.

    Args:
        si_unit: SI unit (e.g., "m", "s")
        quantity: Type of quantity ("length", "time", "temperature")

    Returns:
        Tuple of (display_unit, conversion_factor)
    """
    if quantity == "length":
        # For typical film thicknesses, use µm
        return ("µm", 1e6)
    elif quantity == "time":
        # For typical contact times, use days
        return ("days", 1/86400)
    elif quantity == "temperature":
        return ("°C", 1)  # Already in preferred unit usually
    else:
        return (si_unit, 1.0)


# =============================================================================
# Validation Functions
# =============================================================================

def is_valid_length_unit(unit: str) -> bool:
    """Check if a unit is a valid length unit."""
    return unit.lower().strip() in LENGTH_TO_SI


def is_valid_time_unit(unit: str) -> bool:
    """Check if a unit is a valid time unit."""
    return unit.lower().strip() in TIME_TO_SI


def is_valid_temperature_unit(unit: str) -> bool:
    """Check if a unit is a valid temperature unit."""
    norm = unit.strip().lower().replace("deg", "")
    return norm in ("c", "f", "k", "celsius", "fahrenheit", "kelvin")


def get_available_units() -> Dict[str, list]:
    """Get all available units by category."""
    return {
        "length": list(LENGTH_TO_SI.keys()),
        "time": list(TIME_TO_SI.keys()),
        "temperature": list(TEMPERATURE_UNITS),
        "concentration": CONCENTRATION_UNITS,
    }


# =============================================================================
# ValueWithUnit Helper
# =============================================================================

class UnitValue:
    """
    A numeric value with explicit unit.

    Provides easy conversion and formatting.
    """

    def __init__(self, value: float, unit: str):
        self.value = value
        self.unit = unit

    def to_si(self) -> float:
        """Convert to SI units."""
        return convert_to_si(self.value, self.unit)

    def to_unit(self, target_unit: str) -> float:
        """Convert to a specific unit."""
        unit_lower = self.unit.lower().strip()
        target_lower = target_unit.lower().strip()

        if unit_lower in LENGTH_TO_SI:
            return convert_length(self.value, self.unit, target_unit)
        elif unit_lower in TIME_TO_SI:
            return convert_time(self.value, self.unit, target_unit)
        else:
            return convert_temperature(self.value, self.unit, target_unit)

    def __repr__(self) -> str:
        return f"UnitValue({self.value}, '{self.unit}')"

    def __str__(self) -> str:
        return format_value_with_unit(self.value, self.unit)
