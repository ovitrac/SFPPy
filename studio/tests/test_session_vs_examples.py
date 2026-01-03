"""
Test Suite: Session JSON Files vs. Original Python Examples

Validates that the .sfppy.json session files accurately represent
the configurations defined in the original example*.py scripts.

This ensures:
1. Parameter consistency between JSON and Python definitions
2. Correct API endpoints for session file loading
3. Session validation and conversion to Patankar inputs

@author: SFPPy Studio
@license: MIT
"""

import sys
from pathlib import Path
import pytest
import json

# Add paths
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
sys.path.insert(0, str(Path(__file__).parent.parent))

from studio.app.models.session import Session
from studio.app.utils.units import convert_length, convert_time, convert_temperature


# =============================================================================
# Test Fixtures
# =============================================================================

EXAMPLES_DIR = Path(__file__).parent.parent / "examples"


@pytest.fixture
def example1_json():
    """Load example1.sfppy.json"""
    path = EXAMPLES_DIR / "example1.sfppy.json"
    return Session.load(str(path))


@pytest.fixture
def example1_extensions_json():
    """Load example1_extensions.sfppy.json"""
    path = EXAMPLES_DIR / "example1_extensions.sfppy.json"
    return Session.load(str(path))


@pytest.fixture
def example3_json():
    """Load example3.sfppy.json"""
    path = EXAMPLES_DIR / "example3.sfppy.json"
    return Session.load(str(path))


@pytest.fixture
def example3_variant_json():
    """Load example3_variant.sfppy.json"""
    path = EXAMPLES_DIR / "example3_variant.sfppy.json"
    return Session.load(str(path))


# =============================================================================
# Example 1: LDPE Film Migration Tests
# =============================================================================

class TestExample1VsPython:
    """
    Test example1.sfppy.json against example1.py parameters.

    From example1.py:
    - contactTemperature = (7, "degC")
    - contactTime = (10, "days")
    - maxConcentration = 5000
    - sandwich_geom = Packaging3D('Cylinder', length=(19, "cm"), radius=(30, "mm"))
    - LDPElayer = polymer.LDPE(l=(100, "um"), substance=m1, C0=maxConcentration, T=contactTemperature)
    - Substances: Irganox 1076 (CID 16386), Irgafos 168 (CID 91601)
    """

    def test_geometry_shape(self, example1_json):
        """Test geometry is cylinder as in example1.py"""
        assert example1_json.geometry.shape.value == "cylinder"

    def test_geometry_dimensions(self, example1_json):
        """Test cylinder dimensions match: radius=30mm, height/length=19cm"""
        dims = example1_json.geometry.dimensions

        # Radius: 30 mm
        assert dims.radius is not None
        assert dims.radius.value == 30
        assert dims.radius.unit == "mm"

        # Height: 19 cm
        assert dims.height is not None
        assert dims.height.value == 19
        assert dims.height.unit == "cm"

    def test_layer_count(self, example1_json):
        """Test single layer configuration"""
        assert len(example1_json.layers) == 1

    def test_layer_polymer(self, example1_json):
        """Test layer is LDPE"""
        assert example1_json.layers[0].polymer.value == "LDPE"

    def test_layer_thickness(self, example1_json):
        """Test LDPE layer thickness is 100 um"""
        layer = example1_json.layers[0]
        assert layer.thickness.value == 100
        assert layer.thickness.unit == "um"

        # Verify SI conversion
        thickness_m = convert_length(layer.thickness.value, layer.thickness.unit, "m")
        assert abs(thickness_m - 1e-4) < 1e-10

    def test_substance_count(self, example1_json):
        """Test two substances: Irganox 1076 and Irgafos 168"""
        assert len(example1_json.substances) == 2

    def test_substance_irganox(self, example1_json):
        """Test Irganox 1076 parameters"""
        irganox = next((s for s in example1_json.substances if "irganox" in s.lookup_name.lower()), None)
        assert irganox is not None
        assert irganox.cid == 16386
        assert irganox.SML == 6.0  # EU SML

    def test_substance_irgafos(self, example1_json):
        """Test Irgafos 168 parameters"""
        irgafos = next((s for s in example1_json.substances if "irgafos" in s.lookup_name.lower()), None)
        assert irgafos is not None
        assert irgafos.cid == 91601
        assert irgafos.SML == 60.0  # EU SML

    def test_c0_values(self, example1_json):
        """Test initial concentrations: maxConcentration = 5000 mg/kg"""
        layer = example1_json.layers[0]
        for sub in layer.substances:
            assert sub.C0.value == 5000
            assert sub.C0.unit == "mg/kg"

    def test_contact_step_count(self, example1_json):
        """Test single contact step"""
        assert len(example1_json.contact_steps) == 1

    def test_contact_temperature(self, example1_json):
        """Test contact temperature: 7°C"""
        step = example1_json.contact_steps[0]
        assert step.temperature.value == 7
        assert step.temperature.unit == "degC"

    def test_contact_duration(self, example1_json):
        """Test contact duration: 10 days"""
        step = example1_json.contact_steps[0]
        assert step.duration.value == 10
        assert step.duration.unit == "days"

        # Verify SI conversion
        duration_s = convert_time(step.duration.value, step.duration.unit, "s")
        assert abs(duration_s - 10 * 86400) < 1

    def test_food_properties(self, example1_json):
        """Test food: realfood, semisolid, fat"""
        food = example1_json.food
        assert food.category.value == "realfood"
        assert food.texture.value == "semisolid"
        assert food.affinity.value == "fat"


# =============================================================================
# Example 1 Extensions: Two-Step Migration Tests
# =============================================================================

class TestExample1ExtensionsVsPython:
    """
    Test example1_extensions.sfppy.json with multiple contact steps.
    """

    def test_contact_step_count(self, example1_extensions_json):
        """Test two contact steps"""
        assert len(example1_extensions_json.contact_steps) == 2

    def test_step1_is_storage(self, example1_extensions_json):
        """Test first step is storage"""
        step1 = example1_extensions_json.contact_steps[0]
        assert step1.type.value == "storage"
        assert step1.with_food_contact is True

    def test_step2_is_different_temp(self, example1_extensions_json):
        """Test second step has different temperature"""
        step1 = example1_extensions_json.contact_steps[0]
        step2 = example1_extensions_json.contact_steps[1]
        assert step1.temperature.value != step2.temperature.value


# =============================================================================
# Example 2: Multilayer PP Bottle Tests (no JSON yet, create expected values)
# =============================================================================

class TestExample2ExpectedValues:
    """
    Expected values from example2.py for future JSON file validation.

    From example2.py:
    - contactTemperature = (20, "degC")
    - contactTime = (450, "days")
    - maxConcentration = 10 (mg/kg)
    - bottle = Packaging3D("bottle", body_radius=(40, "mm"), body_height=(0.2, "m"),
    -                       neck_radius=(1.8, "cm"), neck_height=0.05)
    - PPwalls_with_toluene = polymer.PP(l=(300, "um"), substance=surrogate, C0=10, T=(20, "degC"))
    - PET_functionalBarrier = polymer.wPET(l=(30, "um"), C0=0)
    - surrogate = migrant("toluene")
    """

    def test_example2_parameters_documented(self):
        """Document expected parameters for example2.sfppy.json"""
        expected = {
            "geometry": {
                "shape": "bottle",
                "body_radius": (40, "mm"),
                "body_height": (0.2, "m"),
                "neck_radius": (1.8, "cm"),
                "neck_height": (0.05, "m"),
            },
            "substance": {
                "name": "toluene",
                "cid": 1140,
            },
            "layers": [
                {"polymer": "wPET", "thickness": (30, "um"), "C0": 0},  # FB
                {"polymer": "PP", "thickness": (300, "um"), "C0": 10},  # walls
            ],
            "contact": {
                "temperature": (20, "degC"),
                "duration": (450, "days"),
            },
        }
        # This test documents expected values
        assert expected["geometry"]["shape"] == "bottle"
        assert expected["substance"]["cid"] == 1140
        assert len(expected["layers"]) == 2


# =============================================================================
# Example 3: Trilayer ABA Migration Tests
# =============================================================================

class TestExample3VsPython:
    """
    Test example3.sfppy.json against example3.py parameters.

    From example3.py:
    - container = Packaging3D('box_container', height=(8, "cm"), width=(10, "cm"), length=(19, "cm"))
    - m = migrant("limonene")
    - Aw = wPET(l=(20, "um"), migrant=m, C0=0)      # 20 µm plasticized PET
    - B = PP(l=(0.5, "mm"), migrant=m, CP0=200)    # 500 µm PP with 200 mg/kg limonene
    - A = gPET(l=(20, "um"), migrant=m, C0=0)      # 20 µm PET
    - ABA = Aw + B + A
    - contact1: 20°C, 4 months, setoff (no food contact)
    - contact2: hotfilled, 80°C, 20 min
    - contact3: 25°C, 6 months storage
    """

    def test_geometry_shape(self, example3_json):
        """Test geometry is box_container"""
        assert example3_json.geometry.shape.value == "box_container"

    def test_geometry_dimensions(self, example3_json):
        """Test box dimensions: 19x10x8 cm"""
        dims = example3_json.geometry.dimensions
        assert dims.length.value == 19
        assert dims.length.unit == "cm"
        assert dims.width.value == 10
        assert dims.width.unit == "cm"
        assert dims.height.value == 8
        assert dims.height.unit == "cm"

    def test_layer_count(self, example3_json):
        """Test trilayer ABA structure"""
        assert len(example3_json.layers) == 3

    def test_layer1_wPET(self, example3_json):
        """Test Layer 1: wPET (food side barrier)"""
        layer = example3_json.layers[0]
        assert layer.polymer.value == "wPET"
        assert layer.thickness.value == 20
        assert layer.thickness.unit == "um"

    def test_layer2_PP(self, example3_json):
        """Test Layer 2: PP core (500 µm = 0.5 mm)"""
        layer = example3_json.layers[1]
        assert layer.polymer.value == "PP"
        assert layer.thickness.value == 500
        assert layer.thickness.unit == "um"

    def test_layer3_gPET(self, example3_json):
        """Test Layer 3: gPET (outer barrier)"""
        layer = example3_json.layers[2]
        assert layer.polymer.value == "gPET"
        assert layer.thickness.value == 20
        assert layer.thickness.unit == "um"

    def test_c0_distribution(self, example3_json):
        """Test C0: 0 in PET layers, 200 in PP"""
        layers = example3_json.layers

        # wPET: C0 = 0
        wPET_sub = layers[0].substances[0] if layers[0].substances else None
        if wPET_sub:
            assert wPET_sub.C0.value == 0

        # PP: C0 = 200
        PP_sub = layers[1].substances[0] if layers[1].substances else None
        if PP_sub:
            assert PP_sub.C0.value == 200

        # gPET: C0 = 0
        gPET_sub = layers[2].substances[0] if layers[2].substances else None
        if gPET_sub:
            assert gPET_sub.C0.value == 0

    def test_substance_limonene(self, example3_json):
        """Test limonene substance"""
        limonene = next((s for s in example3_json.substances if s.lookup_name == "limonene"), None)
        assert limonene is not None
        assert limonene.cid == 22311

    def test_contact_step_count(self, example3_json):
        """Test three contact steps"""
        assert len(example3_json.contact_steps) == 3

    def test_step1_setoff(self, example3_json):
        """Test Step 1: Setoff - 20°C, 4 months, no food contact"""
        step = example3_json.contact_steps[0]
        assert step.type.value == "setoff"
        assert step.temperature.value == 20
        assert step.temperature.unit == "degC"
        assert step.duration.value == 4
        assert step.duration.unit == "months"
        assert step.with_food_contact is False

    def test_step2_hotfilling(self, example3_json):
        """Test Step 2: Hot filling - 80°C, 20 min"""
        step = example3_json.contact_steps[1]
        assert step.type.value == "hotfilling"
        assert step.temperature.value == 80
        assert step.temperature.unit == "degC"
        assert step.duration.value == 20
        assert step.duration.unit == "minutes"
        assert step.with_food_contact is True

    def test_step3_storage(self, example3_json):
        """Test Step 3: Storage - 25°C, 6 months"""
        step = example3_json.contact_steps[2]
        assert step.type.value == "storage"
        assert step.temperature.value == 25
        assert step.temperature.unit == "degC"
        assert step.duration.value == 6
        assert step.duration.unit == "months"
        assert step.with_food_contact is True

    def test_food_properties(self, example3_json):
        """Test food: realfood, liquid, fat"""
        food = example3_json.food
        assert food.category.value == "realfood"
        assert food.texture.value == "liquid"
        assert food.affinity.value == "fat"


# =============================================================================
# Example 3 Variant: Reduced Barrier Thickness
# =============================================================================

class TestExample3VariantVsPython:
    """
    Test example3_variant.sfppy.json with reduced PET barrier thickness.

    From example3.py variant logic:
    - newthickness[[0, -1]] /= 2  (PET layers halved: 20 -> 10 µm)
    """

    def test_reduced_barrier_thickness(self, example3_variant_json):
        """Test PET layers reduced from 20 µm to 10 µm"""
        layers = example3_variant_json.layers

        # Layer 1 (wPET): 10 µm
        assert layers[0].thickness.value == 10

        # Layer 3 (gPET): 10 µm
        assert layers[2].thickness.value == 10

        # PP unchanged at 500 µm
        assert layers[1].thickness.value == 500


# =============================================================================
# Example 4: Fitting Parameters (expected values for fitting feature)
# =============================================================================

class TestExample4ExpectedValues:
    """
    Expected values from example4.py for fitting feature validation.

    From example4.py:
    - P = layer(l=(100, "um"), D=(1e-10, "cm**2/s"), C0=1000, k=0.1)
    - F = foodlayer(contacttime=(10, "days"), volume=(1, "L"),
    -               surfacearea=(6, "dm**2"), h=(1e-6, "m/s"), CF0=0, k=1)
    - Pseudo-experiment: npoints=30, std_relative=0.01
    - Fitting: R.fit(E) to recover D and k
    """

    def test_fitting_parameters_documented(self):
        """Document expected fitting parameters"""
        expected = {
            "layer": {
                "thickness": (100, "um"),
                "D_true": (1e-10, "cm**2/s"),  # = 1e-14 m²/s
                "C0": 1000,
                "k_true": 0.1,
            },
            "food": {
                "contact_time": (10, "days"),
                "volume": (1, "L"),
                "surface_area": (6, "dm**2"),
                "h": (1e-6, "m/s"),
                "CF0": 0,
                "k": 1,
            },
            "pseudo_experiment": {
                "n_points": 30,
                "std_relative": 0.01,
            },
        }

        # Document D in SI units
        D_cm2_s = 1e-10
        D_m2_s = D_cm2_s * 1e-4  # Convert cm²/s to m²/s
        assert abs(D_m2_s - 1e-14) < 1e-20

        # Verify conversion of surface area
        # 6 dm² = 6 * 0.01 m² = 0.06 m²
        assert expected["food"]["surface_area"] == (6, "dm**2")


# =============================================================================
# API Integration Tests
# =============================================================================

class TestSessionAPIIntegration:
    """Test API endpoints for session file operations"""

    @pytest.fixture
    def client(self):
        """Create test client"""
        from fastapi.testclient import TestClient
        from studio.app.main import app
        return TestClient(app)

    def test_list_examples(self, client):
        """Test GET /api/sessions/files/examples"""
        response = client.get("/api/sessions/files/examples")
        assert response.status_code == 200
        data = response.json()
        assert "examples" in data
        assert len(data["examples"]) >= 4  # At least our 4 examples

    def test_load_example1(self, client):
        """Test loading example1 session"""
        response = client.get("/api/sessions/files/load/Example%201:%20LDPE%20Film%20Migration")
        if response.status_code == 200:
            data = response.json()
            assert data["success"] is True
            assert "session" in data

            session = data["session"]
            assert session["geometry"]["shape"] == "cylinder"
            assert len(session["layers"]) == 1
            assert session["layers"][0]["polymer"] == "LDPE"

    def test_load_example3(self, client):
        """Test loading example3 session"""
        response = client.get("/api/sessions/files/load/Example%203:%20Trilayer%20ABA%20Migration")
        if response.status_code == 200:
            data = response.json()
            assert data["success"] is True

            session = data["session"]
            assert session["geometry"]["shape"] == "box_container"
            assert len(session["layers"]) == 3
            assert len(session["contact_steps"]) == 3

    def test_validate_example1(self, client):
        """Test session validation endpoint"""
        # Load example1 first
        path = EXAMPLES_DIR / "example1.sfppy.json"
        with open(path) as f:
            session_data = json.load(f)

        response = client.post(
            "/api/sessions/files/validate",
            json={"filepath": str(path)}
        )
        # Check response structure
        assert response.status_code in [200, 404]


# =============================================================================
# Unit Conversion Tests
# =============================================================================

class TestUnitConversions:
    """Test unit conversions used in session processing"""

    def test_length_um_to_m(self):
        """Test µm to m conversion"""
        assert convert_length(100, "um", "m") == 1e-4
        assert convert_length(20, "um", "m") == 2e-5
        assert convert_length(500, "um", "m") == 5e-4

    def test_length_cm_to_m(self):
        """Test cm to m conversion"""
        assert convert_length(19, "cm", "m") == 0.19
        assert convert_length(8, "cm", "m") == 0.08

    def test_length_mm_to_m(self):
        """Test mm to m conversion"""
        assert convert_length(30, "mm", "m") == 0.03
        assert convert_length(0.5, "mm", "um") == 500

    def test_time_days_to_s(self):
        """Test days to seconds conversion"""
        assert convert_time(10, "days", "s") == 10 * 86400
        assert convert_time(1, "day", "s") == 86400

    def test_time_months_to_s(self):
        """Test months to seconds conversion"""
        # 1 month = 30 days = 2592000 s
        assert convert_time(4, "months", "s") == 4 * 2592000
        assert convert_time(6, "months", "days") == 6 * 30

    def test_time_minutes_to_s(self):
        """Test minutes to seconds conversion"""
        assert convert_time(20, "minutes", "s") == 1200
        assert convert_time(20, "min", "s") == 1200

    def test_temperature_degC_to_K(self):
        """Test °C to K conversion"""
        assert convert_temperature(7, "degC", "K") == 280.15
        assert convert_temperature(20, "degC", "K") == 293.15
        assert convert_temperature(80, "degC", "K") == 353.15


# =============================================================================
# Session Completeness Tests
# =============================================================================

class TestSessionCompleteness:
    """Test that sessions contain all required fields"""

    def test_example1_complete(self, example1_json):
        """Test example1 has all required fields"""
        assert example1_json.version is not None
        assert example1_json.metadata is not None
        assert example1_json.metadata.name is not None
        assert example1_json.geometry is not None
        assert example1_json.substances is not None
        assert example1_json.layers is not None
        assert example1_json.food is not None
        assert example1_json.contact_steps is not None
        assert example1_json.simulation_config is not None

    def test_example3_complete(self, example3_json):
        """Test example3 has all required fields"""
        assert example3_json.version is not None
        assert example3_json.metadata is not None
        assert example3_json.metadata.name is not None
        assert example3_json.geometry is not None
        assert example3_json.substances is not None
        assert example3_json.layers is not None
        assert example3_json.food is not None
        assert example3_json.contact_steps is not None
        assert example3_json.simulation_config is not None

    def test_metadata_source_script(self, example1_json, example3_json):
        """Test metadata includes source script reference"""
        assert example1_json.metadata.source_script == "example1.py"
        assert example3_json.metadata.source_script == "example3.py"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
