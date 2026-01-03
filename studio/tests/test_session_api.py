"""
API Tests for SFPPy Studio Session Endpoints

Tests the REST API endpoints for session file operations.
Validates each example JSON file can be loaded and processed.

@author: SFPPy Studio
@license: MIT
"""

import sys
from pathlib import Path

# Add paths
sys.path.insert(0, str(Path(__file__).parent.parent))

# Test runner without pytest dependency
import json
import traceback
from typing import Dict, Any, List, Tuple

# Expected values from original example scripts for comparison
EXPECTED_VALUES = {
    "example1": {
        "name": "Example 1: LDPE Film Migration",
        "substances": {
            "count": 2,
            "ids": ["m1", "m2"],
            "names": ["Irganox 1076", "Irgafos 168"],
        },
        "layers": {
            "count": 1,
            "polymers": ["LDPE"],
            "thicknesses_um": [100],
        },
        "contact_steps": {
            "count": 1,
            "temperatures_C": [7],
            "durations_days": [10],
        },
        "geometry": {
            "shape": "cylinder",
            "radius_mm": 30,
            "height_cm": 19,
        },
    },
    "example1_extensions": {
        "name": "Example 1 Extensions: Two-Step Migration",
        "substances": {
            "count": 1,
            "ids": ["m1"],
        },
        "layers": {
            "count": 1,
            "polymers": ["LDPE"],
        },
        "contact_steps": {
            "count": 2,
            "types": ["storage", "ambient"],
            "temperatures_C": [7, 25],
            "durations": [(10, "days"), (4, "hours")],
        },
    },
    "example3": {
        "name": "Example 3: Trilayer ABA Migration",
        "substances": {
            "count": 2,
            "ids": ["limonene", "toluene"],
        },
        "layers": {
            "count": 3,
            "polymers": ["wPET", "PP", "gPET"],
            "thicknesses_um": [20, 500, 20],
        },
        "contact_steps": {
            "count": 3,
            "types": ["setoff", "hotfilling", "storage"],
            "temperatures_C": [20, 80, 25],
        },
        "geometry": {
            "shape": "box_container",
        },
    },
    "example3_variant": {
        "name": "Example 3 Variant: Reduced Barrier Thickness",
        "layers": {
            "count": 3,
            "thicknesses_um": [10, 500, 10],  # Reduced barrier
        },
    },
}


def run_test(test_name: str, test_func) -> Tuple[bool, str]:
    """Run a single test and return result."""
    try:
        test_func()
        return True, ""
    except AssertionError as e:
        return False, str(e)
    except Exception as e:
        return False, f"{type(e).__name__}: {e}\n{traceback.format_exc()}"


def assert_eq(actual, expected, msg=""):
    """Assert equality with detailed message."""
    if actual != expected:
        raise AssertionError(f"{msg}: expected {expected}, got {actual}")


def assert_in(item, collection, msg=""):
    """Assert item is in collection."""
    if item not in collection:
        raise AssertionError(f"{msg}: {item} not found in {collection}")


# =============================================================================
# Test: Load Session Models
# =============================================================================

def test_load_session_models():
    """Test that session models can be imported."""
    from app.models.session import (
        Session, Metadata, Geometry, Substance, Layer, ContactStep, Food
    )
    # Just check imports work
    assert Session is not None


def test_load_utils():
    """Test that utilities can be imported."""
    from app.utils.units import convert_length, convert_time, convert_temperature
    from app.utils.session_io import load_session_file, validate_session


# =============================================================================
# Test: Example 1 - Multiple Substances
# =============================================================================

def test_example1_load():
    """Test loading example1.sfppy.json."""
    from app.utils.session_io import load_session_file
    examples_dir = Path(__file__).parent.parent / "examples"
    session = load_session_file(str(examples_dir / "example1.sfppy.json"))

    expected = EXPECTED_VALUES["example1"]
    assert_eq(session.metadata.name, expected["name"], "Session name")


def test_example1_substances():
    """Test example1 has correct substances."""
    from app.utils.session_io import load_session_file
    examples_dir = Path(__file__).parent.parent / "examples"
    session = load_session_file(str(examples_dir / "example1.sfppy.json"))

    expected = EXPECTED_VALUES["example1"]["substances"]
    assert_eq(len(session.substances), expected["count"], "Substance count")

    substance_ids = [s.id for s in session.substances]
    for sid in expected["ids"]:
        assert_in(sid, substance_ids, "Substance ID")

    # Check substance names
    for sub in session.substances:
        if sub.id == "m1":
            assert_eq(sub.properties.name, "Irganox 1076", "m1 name")
            assert_eq(sub.properties.mw, 530.9, "m1 MW")
            assert_eq(sub.SML, 6.0, "m1 SML")
        elif sub.id == "m2":
            assert_eq(sub.properties.name, "Irgafos 168", "m2 name")
            assert_eq(sub.properties.mw, 646.9, "m2 MW")
            assert_eq(sub.SML, 60.0, "m2 SML")


def test_example1_layers():
    """Test example1 has correct layer structure."""
    from app.utils.session_io import load_session_file
    examples_dir = Path(__file__).parent.parent / "examples"
    session = load_session_file(str(examples_dir / "example1.sfppy.json"))

    expected = EXPECTED_VALUES["example1"]["layers"]
    assert_eq(len(session.layers), expected["count"], "Layer count")
    assert_eq(session.layers[0].polymer.value, expected["polymers"][0], "Layer polymer")
    assert_eq(session.layers[0].thickness.value, expected["thicknesses_um"][0], "Layer thickness")
    assert_eq(session.layers[0].thickness.unit, "um", "Layer thickness unit")


def test_example1_contact_step():
    """Test example1 has correct contact step."""
    from app.utils.session_io import load_session_file
    examples_dir = Path(__file__).parent.parent / "examples"
    session = load_session_file(str(examples_dir / "example1.sfppy.json"))

    expected = EXPECTED_VALUES["example1"]["contact_steps"]
    assert_eq(len(session.contact_steps), expected["count"], "Step count")
    assert_eq(session.contact_steps[0].temperature.value, expected["temperatures_C"][0], "Temperature")
    assert_eq(session.contact_steps[0].duration.value, expected["durations_days"][0], "Duration")


def test_example1_geometry():
    """Test example1 has correct geometry."""
    from app.utils.session_io import load_session_file
    examples_dir = Path(__file__).parent.parent / "examples"
    session = load_session_file(str(examples_dir / "example1.sfppy.json"))

    expected = EXPECTED_VALUES["example1"]["geometry"]
    assert_eq(session.geometry.shape.value, expected["shape"], "Shape")
    assert_eq(session.geometry.dimensions.radius.value, expected["radius_mm"], "Radius")
    assert_eq(session.geometry.dimensions.height.value, expected["height_cm"], "Height")


def test_example1_validation():
    """Test example1 passes validation."""
    from app.utils.session_io import load_session_file, validate_session
    examples_dir = Path(__file__).parent.parent / "examples"
    session = load_session_file(str(examples_dir / "example1.sfppy.json"))

    is_valid, errors, warnings = validate_session(session)
    assert is_valid, f"Validation failed: {errors}"


# =============================================================================
# Test: Example 1 Extensions - Multiple Contact Steps
# =============================================================================

def test_example1_ext_load():
    """Test loading example1_extensions.sfppy.json."""
    from app.utils.session_io import load_session_file
    examples_dir = Path(__file__).parent.parent / "examples"
    session = load_session_file(str(examples_dir / "example1_extensions.sfppy.json"))

    expected = EXPECTED_VALUES["example1_extensions"]
    assert_eq(session.metadata.name, expected["name"], "Session name")


def test_example1_ext_multiple_steps():
    """Test example1_extensions has two contact steps."""
    from app.utils.session_io import load_session_file
    examples_dir = Path(__file__).parent.parent / "examples"
    session = load_session_file(str(examples_dir / "example1_extensions.sfppy.json"))

    expected = EXPECTED_VALUES["example1_extensions"]["contact_steps"]
    assert_eq(len(session.contact_steps), expected["count"], "Step count")

    # Check step types
    for i, step in enumerate(session.contact_steps):
        assert_eq(step.type.value, expected["types"][i], f"Step {i+1} type")
        assert_eq(step.temperature.value, expected["temperatures_C"][i], f"Step {i+1} temperature")


def test_example1_ext_step_durations():
    """Test example1_extensions has correct step durations with different units."""
    from app.utils.session_io import load_session_file
    examples_dir = Path(__file__).parent.parent / "examples"
    session = load_session_file(str(examples_dir / "example1_extensions.sfppy.json"))

    # Step 1: 10 days
    assert_eq(session.contact_steps[0].duration.value, 10, "Step 1 duration value")
    assert_eq(session.contact_steps[0].duration.unit, "days", "Step 1 duration unit")

    # Step 2: 4 hours
    assert_eq(session.contact_steps[1].duration.value, 4, "Step 2 duration value")
    assert_eq(session.contact_steps[1].duration.unit, "hours", "Step 2 duration unit")


# =============================================================================
# Test: Example 3 - Multiple Layers (ABA Structure)
# =============================================================================

def test_example3_load():
    """Test loading example3.sfppy.json."""
    from app.utils.session_io import load_session_file
    examples_dir = Path(__file__).parent.parent / "examples"
    session = load_session_file(str(examples_dir / "example3.sfppy.json"))

    expected = EXPECTED_VALUES["example3"]
    assert_eq(session.metadata.name, expected["name"], "Session name")


def test_example3_aba_structure():
    """Test example3 has correct ABA trilayer structure."""
    from app.utils.session_io import load_session_file
    examples_dir = Path(__file__).parent.parent / "examples"
    session = load_session_file(str(examples_dir / "example3.sfppy.json"))

    expected = EXPECTED_VALUES["example3"]["layers"]
    assert_eq(len(session.layers), expected["count"], "Layer count")

    # Check ABA structure: wPET / PP / gPET
    for i, layer in enumerate(session.layers):
        assert_eq(layer.polymer.value, expected["polymers"][i], f"Layer {i+1} polymer")
        assert_eq(layer.thickness.value, expected["thicknesses_um"][i], f"Layer {i+1} thickness")


def test_example3_three_contact_steps():
    """Test example3 has three contact steps with correct types."""
    from app.utils.session_io import load_session_file
    examples_dir = Path(__file__).parent.parent / "examples"
    session = load_session_file(str(examples_dir / "example3.sfppy.json"))

    expected = EXPECTED_VALUES["example3"]["contact_steps"]
    assert_eq(len(session.contact_steps), expected["count"], "Step count")

    # Check step types: setoff, hotfilling, storage
    for i, step in enumerate(session.contact_steps):
        assert_eq(step.type.value, expected["types"][i], f"Step {i+1} type")
        assert_eq(step.temperature.value, expected["temperatures_C"][i], f"Step {i+1} temperature")


def test_example3_setoff_no_food_contact():
    """Test example3 setoff step has no food contact."""
    from app.utils.session_io import load_session_file
    examples_dir = Path(__file__).parent.parent / "examples"
    session = load_session_file(str(examples_dir / "example3.sfppy.json"))

    # Setoff step should have with_food_contact = False
    setoff_step = session.contact_steps[0]
    assert_eq(setoff_step.type.value, "setoff", "First step is setoff")
    assert_eq(setoff_step.with_food_contact, False, "Setoff has no food contact")


def test_example3_box_geometry():
    """Test example3 has box_container geometry."""
    from app.utils.session_io import load_session_file
    examples_dir = Path(__file__).parent.parent / "examples"
    session = load_session_file(str(examples_dir / "example3.sfppy.json"))

    assert_eq(session.geometry.shape.value, "box_container", "Geometry shape")


# =============================================================================
# Test: Example 3 Variant - Reduced Barrier Thickness
# =============================================================================

def test_example3_variant_load():
    """Test loading example3_variant.sfppy.json."""
    from app.utils.session_io import load_session_file
    examples_dir = Path(__file__).parent.parent / "examples"
    session = load_session_file(str(examples_dir / "example3_variant.sfppy.json"))

    # Check it loaded
    assert session is not None


def test_example3_variant_reduced_thickness():
    """Test example3_variant has reduced PET barrier thickness."""
    from app.utils.session_io import load_session_file
    examples_dir = Path(__file__).parent.parent / "examples"
    session = load_session_file(str(examples_dir / "example3_variant.sfppy.json"))

    expected = EXPECTED_VALUES["example3_variant"]["layers"]

    # PET layers should be 10 µm instead of 20 µm
    assert_eq(session.layers[0].thickness.value, expected["thicknesses_um"][0], "Layer 1 thickness (reduced)")
    assert_eq(session.layers[1].thickness.value, expected["thicknesses_um"][1], "Layer 2 thickness (PP same)")
    assert_eq(session.layers[2].thickness.value, expected["thicknesses_um"][2], "Layer 3 thickness (reduced)")


# =============================================================================
# Test: Unit Conversions
# =============================================================================

def test_unit_conversion_length():
    """Test length unit conversions."""
    from app.utils.units import convert_length

    # um to m
    assert abs(convert_length(100, "um", "m") - 1e-4) < 1e-10

    # mm to um
    assert abs(convert_length(1, "mm", "um") - 1000) < 1e-6

    # cm to mm
    assert abs(convert_length(10, "cm", "mm") - 100) < 1e-6


def test_unit_conversion_time():
    """Test time unit conversions."""
    from app.utils.units import convert_time

    # days to seconds
    assert abs(convert_time(1, "days", "s") - 86400) < 1

    # hours to minutes
    assert abs(convert_time(1, "hours", "min") - 60) < 0.1

    # months to days
    assert abs(convert_time(1, "months", "days") - 30) < 1


def test_unit_conversion_temperature():
    """Test temperature unit conversions."""
    from app.utils.units import convert_temperature

    # degC to K
    assert abs(convert_temperature(25, "degC", "K") - 298.15) < 0.01

    # K to degC
    assert abs(convert_temperature(273.15, "K", "degC") - 0) < 0.01


# =============================================================================
# Test: Session Validation
# =============================================================================

def test_validation_all_examples():
    """Test all examples pass validation."""
    from app.utils.session_io import load_session_file, validate_session
    examples_dir = Path(__file__).parent.parent / "examples"

    for example_name in ["example1", "example1_extensions", "example3", "example3_variant"]:
        session = load_session_file(str(examples_dir / f"{example_name}.sfppy.json"))
        is_valid, errors, warnings = validate_session(session)
        assert is_valid, f"{example_name} validation failed: {errors}"


# =============================================================================
# Test: Serialization Roundtrip
# =============================================================================

def test_roundtrip_all_examples():
    """Test all examples survive JSON roundtrip."""
    from app.utils.session_io import load_session_file
    from app.models.session import Session
    examples_dir = Path(__file__).parent.parent / "examples"

    for example_name in ["example1", "example1_extensions", "example3", "example3_variant"]:
        session = load_session_file(str(examples_dir / f"{example_name}.sfppy.json"))

        # Serialize and deserialize
        json_str = session.to_json()
        reloaded = Session.from_json(json_str)

        assert_eq(reloaded.metadata.name, session.metadata.name, f"{example_name} name")
        assert_eq(len(reloaded.layers), len(session.layers), f"{example_name} layers")
        assert_eq(len(reloaded.contact_steps), len(session.contact_steps), f"{example_name} steps")


# =============================================================================
# Test Runner
# =============================================================================

def run_all_tests() -> Tuple[int, int, List[str]]:
    """Run all tests and return results."""
    tests = [
        # Model imports
        ("test_load_session_models", test_load_session_models),
        ("test_load_utils", test_load_utils),

        # Example 1
        ("test_example1_load", test_example1_load),
        ("test_example1_substances", test_example1_substances),
        ("test_example1_layers", test_example1_layers),
        ("test_example1_contact_step", test_example1_contact_step),
        ("test_example1_geometry", test_example1_geometry),
        ("test_example1_validation", test_example1_validation),

        # Example 1 Extensions
        ("test_example1_ext_load", test_example1_ext_load),
        ("test_example1_ext_multiple_steps", test_example1_ext_multiple_steps),
        ("test_example1_ext_step_durations", test_example1_ext_step_durations),

        # Example 3
        ("test_example3_load", test_example3_load),
        ("test_example3_aba_structure", test_example3_aba_structure),
        ("test_example3_three_contact_steps", test_example3_three_contact_steps),
        ("test_example3_setoff_no_food_contact", test_example3_setoff_no_food_contact),
        ("test_example3_box_geometry", test_example3_box_geometry),

        # Example 3 Variant
        ("test_example3_variant_load", test_example3_variant_load),
        ("test_example3_variant_reduced_thickness", test_example3_variant_reduced_thickness),

        # Unit conversions
        ("test_unit_conversion_length", test_unit_conversion_length),
        ("test_unit_conversion_time", test_unit_conversion_time),
        ("test_unit_conversion_temperature", test_unit_conversion_temperature),

        # Validation
        ("test_validation_all_examples", test_validation_all_examples),

        # Roundtrip
        ("test_roundtrip_all_examples", test_roundtrip_all_examples),
    ]

    passed = 0
    failed = 0
    failures = []

    print("=" * 60)
    print("SFPPy Studio API Tests")
    print("=" * 60)
    print()

    for name, func in tests:
        success, error = run_test(name, func)
        if success:
            print(f"  [PASS] {name}")
            passed += 1
        else:
            print(f"  [FAIL] {name}")
            print(f"         {error.split(chr(10))[0]}")
            failed += 1
            failures.append((name, error))

    print()
    print("=" * 60)
    print(f"SUMMARY: {passed} passed, {failed} failed")
    print("=" * 60)

    if failures:
        print()
        print("FAILURES:")
        for name, error in failures:
            print(f"\n{name}:")
            print(error)

    return passed, failed, failures


if __name__ == "__main__":
    import os
    os.chdir(Path(__file__).parent.parent)
    passed, failed, _ = run_all_tests()
    sys.exit(1 if failed > 0 else 0)
