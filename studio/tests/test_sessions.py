"""
SFPPy Studio Session API Tests

Comprehensive test suite for session-based substance management,
migrantToxtree storage, and concentration unit conversions.

Run with: pytest studio/tests/test_sessions.py -v
"""

import pytest
import sys
from pathlib import Path

# Add paths for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from fastapi.testclient import TestClient
from studio.app.main import app
from studio.app.session import (
    convert_concentration,
    get_concentration_for_computation,
    CONCENTRATION_UNITS,
    TTC_VALUES,
    CFTTC_VALUES,
)

client = TestClient(app)


# ========== UNIT CONVERSION TESTS ==========

class TestConcentrationConversions:
    """Test concentration unit conversion functions."""

    def test_mg_kg_to_kg_m3(self):
        """Test mg/kg to kg/m³ conversion."""
        # 1000 mg/kg at density 920 kg/m³ (LDPE)
        # = 1000 * 1e-6 kg/kg * 920 kg/m³ = 0.92 kg/m³
        result = convert_concentration(1000, "mg/kg", "kg/m3", 920)
        assert result == pytest.approx(0.92, rel=0.01)

    def test_kg_m3_to_mg_kg(self):
        """Test kg/m³ to mg/kg conversion."""
        # 0.92 kg/m³ at density 920 kg/m³
        # = 0.92 / 920 / 1e-6 = 1000 mg/kg
        result = convert_concentration(0.92, "kg/m3", "mg/kg", 920)
        assert result == pytest.approx(1000, rel=0.01)

    def test_ppm_equals_mg_kg(self):
        """Test ppm is equivalent to mg/kg."""
        value = 500
        result_ppm = convert_concentration(value, "ppm", "kg/m3", 1000)
        result_mg_kg = convert_concentration(value, "mg/kg", "kg/m3", 1000)
        assert result_ppm == result_mg_kg

    def test_ppb_equals_ug_kg(self):
        """Test ppb is equivalent to µg/kg."""
        value = 100
        result_ppb = convert_concentration(value, "ppb", "kg/m3", 1000)
        result_ug_kg = convert_concentration(value, "µg/kg", "kg/m3", 1000)
        assert result_ppb == result_ug_kg

    def test_g_kg_to_mg_kg(self):
        """Test g/kg to mg/kg conversion."""
        # 1 g/kg = 1000 mg/kg
        result = convert_concentration(1, "g/kg", "mg/kg", 1000)
        assert result == pytest.approx(1000, rel=0.01)

    def test_kg_kg_is_identity(self):
        """Test kg/kg is the base unit (mass fraction)."""
        # 0.001 kg/kg = 1 g/kg = 1000 mg/kg
        result = convert_concentration(0.001, "kg/kg", "mg/kg", 1000)
        assert result == pytest.approx(1000, rel=0.01)

    def test_ng_kg_to_ppb(self):
        """Test ng/kg to ppb conversion."""
        # 1000 ng/kg = 1 ppb (= 1 µg/kg)
        result = convert_concentration(1000, "ng/kg", "ppb", 1000)
        assert result == pytest.approx(1, rel=0.01)

    def test_round_trip_conversion(self):
        """Test round-trip conversion preserves value."""
        original = 5000
        to_kg_m3 = convert_concentration(original, "mg/kg", "kg/m3", 1000)
        back = convert_concentration(to_kg_m3, "kg/m3", "mg/kg", 1000)
        assert back == pytest.approx(original, rel=0.001)

    def test_same_unit_no_conversion(self):
        """Test same unit returns same value."""
        result = convert_concentration(100, "mg/kg", "mg/kg", 1000)
        assert result == 100

    def test_invalid_unit_raises(self):
        """Test invalid unit raises ValueError."""
        with pytest.raises(ValueError):
            convert_concentration(100, "invalid_unit", "mg/kg", 1000)

    def test_get_concentration_for_computation(self):
        """Test convenience function for computation conversion."""
        result = get_concentration_for_computation(1000, "mg/kg", 920)
        assert result == pytest.approx(0.92, rel=0.01)


class TestConcentrationUnitsAPI:
    """Test concentration units API endpoints."""

    def test_list_concentration_units(self):
        """Test listing available concentration units."""
        response = client.get("/api/sessions/units/concentration")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["default"] == "mg/kg"
        assert len(data["units"]) >= 7
        # Check required units
        codes = [u["code"] for u in data["units"]]
        assert "mg/kg" in codes
        assert "ppm" in codes
        assert "ppb" in codes
        assert "g/kg" in codes
        assert "kg/m3" in codes

    def test_convert_concentration_api(self):
        """Test concentration conversion endpoint."""
        response = client.post(
            "/api/sessions/units/convert",
            params={
                "value": 1000,
                "from_unit": "mg/kg",
                "to_unit": "g/kg",
                "polymer_density": 1000
            }
        )
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["output"]["value"] == pytest.approx(1.0, rel=0.01)
        assert data["output"]["unit"] == "g/kg"

    def test_convert_invalid_unit_api(self):
        """Test conversion with invalid unit."""
        response = client.post(
            "/api/sessions/units/convert",
            params={
                "value": 100,
                "from_unit": "invalid",
                "to_unit": "mg/kg",
                "polymer_density": 1000
            }
        )
        assert response.status_code == 400


# ========== TTC REFERENCE TESTS ==========

class TestTTCReference:
    """Test TTC reference data."""

    def test_ttc_values_defined(self):
        """Test TTC values are properly defined."""
        assert len(TTC_VALUES) == 4
        # Class I (index 1) should be 1.5 µg/kg bw/day
        assert TTC_VALUES[1] == pytest.approx(1.5, rel=0.01)
        # Class III (index 3) should be 30 µg/kg bw/day
        assert TTC_VALUES[3] == pytest.approx(30, rel=0.01)

    def test_cfttc_values_calculated(self):
        """Test CF_TTC values are calculated from TTC."""
        # CF_TTC = TTC * 60 * 1e-3 (60 kg bw, 1 kg food, mg conversion)
        assert len(CFTTC_VALUES) == 4
        for i in range(4):
            expected = TTC_VALUES[i] * 60 * 1e-3
            assert CFTTC_VALUES[i] == pytest.approx(expected, rel=0.01)

    def test_ttc_reference_api(self):
        """Test TTC reference endpoint."""
        response = client.get("/api/sessions/reference/ttc")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "ttc" in data
        assert "cf_ttc" in data
        assert len(data["ttc"]["values"]) == 4
        assert data["ttc"]["unit"] == "µg/kg bw/day"
        assert data["cf_ttc"]["unit"] == "mg/kg food"


# ========== SESSION MANAGEMENT TESTS ==========

class TestSessionManagement:
    """Test session creation, listing, and deletion."""

    def test_create_session(self):
        """Test creating a new session."""
        response = client.post(
            "/api/sessions/create",
            json={"name": "Test Session"}
        )
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "session" in data
        assert "id" in data["session"]
        assert data["session"]["name"] == "Test Session"

    def test_create_session_default_name(self):
        """Test creating session with default name."""
        response = client.post("/api/sessions/create", json={})
        assert response.status_code == 200
        data = response.json()
        assert data["session"]["name"] == "Untitled Simulation"

    def test_list_sessions(self):
        """Test listing all sessions."""
        # Create a session first
        client.post("/api/sessions/create", json={"name": "List Test"})

        response = client.get("/api/sessions/list")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "sessions" in data
        assert data["count"] >= 1

    def test_get_session(self):
        """Test getting session details."""
        # Create session
        create_resp = client.post("/api/sessions/create", json={"name": "Get Test"})
        session_id = create_resp.json()["session"]["id"]

        # Get session
        response = client.get(f"/api/sessions/{session_id}")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["session"]["id"] == session_id
        assert data["session"]["name"] == "Get Test"
        assert "substances" in data["session"]
        assert "layers" in data["session"]
        assert "contact_steps" in data["session"]

    def test_get_nonexistent_session(self):
        """Test getting non-existent session."""
        response = client.get("/api/sessions/nonexistent-id-12345")
        assert response.status_code == 404

    def test_delete_session(self):
        """Test deleting a session."""
        # Create session
        create_resp = client.post("/api/sessions/create", json={"name": "Delete Test"})
        session_id = create_resp.json()["session"]["id"]

        # Delete session
        response = client.delete(f"/api/sessions/{session_id}")
        assert response.status_code == 200

        # Verify deleted
        response = client.get(f"/api/sessions/{session_id}")
        assert response.status_code == 404


# ========== SUBSTANCE MANAGEMENT TESTS ==========

class TestSubstanceManagement:
    """Test substance addition and removal in sessions."""

    @pytest.fixture
    def session_id(self):
        """Create a session for testing."""
        response = client.post("/api/sessions/create", json={"name": "Substance Test"})
        return response.json()["session"]["id"]

    def test_add_substance_by_name(self, session_id):
        """Test adding substance by name."""
        response = client.post(
            f"/api/sessions/{session_id}/substances",
            json={"query": "anisole"}
        )
        # May fail if PubChem not available, skip in that case
        if response.status_code == 404:
            pytest.skip("PubChem lookup not available")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["added"] is True
        assert "substance" in data

    def test_add_substance_by_cas(self, session_id):
        """Test adding substance by CAS number."""
        response = client.post(
            f"/api/sessions/{session_id}/substances",
            json={"query": "100-66-3"}  # Anisole CAS
        )
        if response.status_code == 404:
            pytest.skip("PubChem lookup not available")
        assert response.status_code == 200

    def test_add_unknown_substance(self, session_id):
        """Test adding unknown substance with defaults."""
        response = client.post(
            f"/api/sessions/{session_id}/substances/unknown",
            json={"name": "mystery_compound", "mw": 300, "logP": 4.5}
        )
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["substance"]["D"] == 1e-12  # Default D
        assert data["substance"]["k"] == 1.0    # Default k

    def test_add_duplicate_substance(self, session_id):
        """Test adding same substance twice."""
        # Add unknown substance
        client.post(
            f"/api/sessions/{session_id}/substances/unknown",
            json={"name": "test_dup", "mw": 200}
        )

        # Add again
        response = client.post(
            f"/api/sessions/{session_id}/substances/unknown",
            json={"name": "test_dup", "mw": 200}
        )
        assert response.status_code == 200
        data = response.json()
        assert data["added"] is False
        assert "already in session" in data["message"]

    def test_list_session_substances(self, session_id):
        """Test listing substances in session."""
        # Add a substance
        client.post(
            f"/api/sessions/{session_id}/substances/unknown",
            json={"name": "list_test", "mw": 250}
        )

        response = client.get(f"/api/sessions/{session_id}/substances")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["count"] >= 1

    def test_remove_substance(self, session_id):
        """Test removing substance from session."""
        # Add a substance
        add_resp = client.post(
            f"/api/sessions/{session_id}/substances/unknown",
            json={"name": "remove_test", "mw": 300}
        )
        substance_id = add_resp.json()["substance"]["id"]

        # Remove it
        response = client.delete(f"/api/sessions/{session_id}/substances/{substance_id}")
        assert response.status_code == 200

        # Verify removed
        response = client.get(f"/api/sessions/{session_id}/substances")
        data = response.json()
        substance_ids = [s["id"] for s in data["substances"]]
        assert substance_id not in substance_ids

    def test_get_substance_details(self, session_id):
        """Test getting substance details."""
        # Add a substance
        add_resp = client.post(
            f"/api/sessions/{session_id}/substances/unknown",
            json={"name": "detail_test", "mw": 400, "logP": 6}
        )
        substance_id = add_resp.json()["substance"]["id"]

        response = client.get(f"/api/sessions/{session_id}/substances/{substance_id}")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["substance"]["mw"] == 400
        assert data["substance"]["logP"] == 6


# ========== LAYER MANAGEMENT TESTS ==========

class TestLayerManagement:
    """Test layer addition and removal in sessions."""

    @pytest.fixture
    def session_id(self):
        """Create a session for testing."""
        response = client.post("/api/sessions/create", json={"name": "Layer Test"})
        return response.json()["session"]["id"]

    def test_add_layer(self, session_id):
        """Test adding a layer."""
        response = client.post(
            f"/api/sessions/{session_id}/layers",
            json={"polymer": "LDPE", "thickness": 100, "thickness_unit": "um"}
        )
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["layer"]["index"] == 1
        assert data["layer"]["polymer"] == "LDPE"
        assert data["layer"]["thickness"] == (100, "um")

    def test_add_multiple_layers(self, session_id):
        """Test adding multiple layers with auto-indexing."""
        # Add first layer
        client.post(
            f"/api/sessions/{session_id}/layers",
            json={"polymer": "gPET", "thickness": 20, "thickness_unit": "um"}
        )

        # Add second layer
        response = client.post(
            f"/api/sessions/{session_id}/layers",
            json={"polymer": "LDPE", "thickness": 100, "thickness_unit": "um"}
        )
        assert response.status_code == 200
        data = response.json()
        assert data["layer"]["index"] == 2

    def test_list_layers(self, session_id):
        """Test listing layers."""
        client.post(
            f"/api/sessions/{session_id}/layers",
            json={"polymer": "PP", "thickness": 200, "thickness_unit": "um"}
        )

        response = client.get(f"/api/sessions/{session_id}/layers")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["count"] >= 1

    def test_remove_layer_reindex(self, session_id):
        """Test removing layer re-indexes remaining layers."""
        # Add 3 layers
        for polymer in ["LDPE", "PP", "PET"]:
            client.post(
                f"/api/sessions/{session_id}/layers",
                json={"polymer": polymer, "thickness": 50, "thickness_unit": "um"}
            )

        # Remove middle layer (index 2)
        response = client.delete(f"/api/sessions/{session_id}/layers/2")
        assert response.status_code == 200

        # Check re-indexing
        response = client.get(f"/api/sessions/{session_id}/layers")
        data = response.json()
        assert data["count"] == 2
        indices = [l["index"] for l in data["layers"]]
        assert indices == [1, 2]  # Re-indexed


# ========== CONTACT STEP TESTS ==========

class TestContactStepManagement:
    """Test contact step management in sessions."""

    @pytest.fixture
    def session_id(self):
        """Create a session for testing."""
        response = client.post("/api/sessions/create", json={"name": "Step Test"})
        return response.json()["session"]["id"]

    def test_add_contact_step(self, session_id):
        """Test adding a contact step."""
        response = client.post(
            f"/api/sessions/{session_id}/steps",
            json={
                "temperature": 40,
                "temperature_unit": "degC",
                "duration": 10,
                "duration_unit": "days"
            }
        )
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["step"]["index"] == 1
        assert data["step"]["temperature"] == (40, "degC")
        assert data["step"]["duration"] == (10, "days")
        assert data["step"]["with_food"] is True

    def test_add_setoff_step(self, session_id):
        """Test adding set-off step (no food contact)."""
        response = client.post(
            f"/api/sessions/{session_id}/steps",
            json={
                "temperature": 25,
                "temperature_unit": "degC",
                "duration": 90,
                "duration_unit": "days",
                "with_food": False
            }
        )
        assert response.status_code == 200
        data = response.json()
        assert data["step"]["with_food"] is False

    def test_list_contact_steps(self, session_id):
        """Test listing contact steps."""
        client.post(
            f"/api/sessions/{session_id}/steps",
            json={"temperature": 40, "duration": 10, "duration_unit": "days"}
        )

        response = client.get(f"/api/sessions/{session_id}/steps")
        assert response.status_code == 200
        data = response.json()
        assert data["count"] >= 1


# ========== ASSIGNMENT TESTS ==========

class TestSubstanceAssignment:
    """Test substance-to-layer assignment."""

    @pytest.fixture
    def session_with_data(self):
        """Create session with substance and layer."""
        # Create session
        resp = client.post("/api/sessions/create", json={"name": "Assignment Test"})
        session_id = resp.json()["session"]["id"]

        # Add substance
        resp = client.post(
            f"/api/sessions/{session_id}/substances/unknown",
            json={"name": "test_sub", "mw": 400, "logP": 5}
        )
        substance_id = resp.json()["substance"]["id"]

        # Add layer
        client.post(
            f"/api/sessions/{session_id}/layers",
            json={"polymer": "LDPE", "thickness": 100, "thickness_unit": "um"}
        )

        return session_id, substance_id

    def test_assign_substance_to_layer(self, session_with_data):
        """Test assigning substance to layer with C0."""
        session_id, substance_id = session_with_data

        response = client.post(
            f"/api/sessions/{session_id}/assignments",
            json={
                "substance_id": substance_id,
                "layer_index": 1,
                "C0": 1000,
                "C0_unit": "mg/kg"
            }
        )
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["assignment"]["C0"] == 1000
        assert data["assignment"]["C0_unit"] == "mg/kg"

    def test_assign_with_override(self, session_with_data):
        """Test assigning with D/k overrides."""
        session_id, substance_id = session_with_data

        response = client.post(
            f"/api/sessions/{session_id}/assignments",
            json={
                "substance_id": substance_id,
                "layer_index": 1,
                "C0": 500,
                "C0_unit": "mg/kg",
                "D_override": 5e-13,
                "k_override": 2.0
            }
        )
        assert response.status_code == 200
        data = response.json()
        assert data["assignment"]["D_override"] == 5e-13
        assert data["assignment"]["k_override"] == 2.0

    def test_list_assignments(self, session_with_data):
        """Test listing assignments."""
        session_id, substance_id = session_with_data

        client.post(
            f"/api/sessions/{session_id}/assignments",
            json={"substance_id": substance_id, "layer_index": 1, "C0": 1000}
        )

        response = client.get(f"/api/sessions/{session_id}/assignments")
        assert response.status_code == 200
        data = response.json()
        assert data["count"] >= 1

    def test_remove_assignment(self, session_with_data):
        """Test removing an assignment."""
        session_id, substance_id = session_with_data

        client.post(
            f"/api/sessions/{session_id}/assignments",
            json={"substance_id": substance_id, "layer_index": 1, "C0": 1000}
        )

        response = client.delete(f"/api/sessions/{session_id}/assignments/{substance_id}/1")
        assert response.status_code == 200


# ========== VALIDATION TESTS ==========

class TestSessionValidation:
    """Test session validation for simulation readiness."""

    @pytest.fixture
    def session_id(self):
        """Create a session for testing."""
        response = client.post("/api/sessions/create", json={"name": "Validation Test"})
        return response.json()["session"]["id"]

    def test_empty_session_invalid(self, session_id):
        """Test empty session is not valid for simulation."""
        response = client.get(f"/api/sessions/{session_id}/validate")
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is False
        assert "At least 1 layer" in str(data["errors"])
        assert "At least 1 contact step" in str(data["errors"])
        assert "At least 1 substance" in str(data["errors"])

    def test_complete_session_valid(self, session_id):
        """Test complete session is valid for simulation."""
        # Add substance
        resp = client.post(
            f"/api/sessions/{session_id}/substances/unknown",
            json={"name": "test", "mw": 300}
        )
        substance_id = resp.json()["substance"]["id"]

        # Add layer
        client.post(
            f"/api/sessions/{session_id}/layers",
            json={"polymer": "LDPE", "thickness": 100, "thickness_unit": "um"}
        )

        # Add step
        client.post(
            f"/api/sessions/{session_id}/steps",
            json={"temperature": 40, "duration": 10, "duration_unit": "days"}
        )

        # Add assignment
        client.post(
            f"/api/sessions/{session_id}/assignments",
            json={"substance_id": substance_id, "layer_index": 1, "C0": 1000}
        )

        # Validate
        response = client.get(f"/api/sessions/{session_id}/validate")
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is True
        assert len(data["errors"]) == 0


# ========== TOXTREE DATA TESTS ==========

class TestToxtreeData:
    """Test ToxTree data availability in substances."""

    @pytest.fixture
    def session_id(self):
        """Create a session for testing."""
        response = client.post("/api/sessions/create", json={"name": "ToxTree Test"})
        return response.json()["session"]["id"]

    def test_unknown_substance_no_toxtree(self, session_id):
        """Test unknown substance has no ToxTree data."""
        response = client.post(
            f"/api/sessions/{session_id}/substances/unknown",
            json={"name": "generic", "mw": 500}
        )
        assert response.status_code == 200
        data = response.json()
        substance = data["substance"]
        assert substance["toxtree"]["cramer_class"] is None
        assert substance["toxtree"]["ttc"] is None

    def test_substance_data_structure(self, session_id):
        """Test substance data has expected structure."""
        response = client.post(
            f"/api/sessions/{session_id}/substances/unknown",
            json={"name": "structure_test", "mw": 400, "logP": 5}
        )
        assert response.status_code == 200
        substance = response.json()["substance"]

        # Check all expected fields
        assert "id" in substance
        assert "name" in substance
        assert "mw" in substance
        assert "logP" in substance
        assert "regulatory" in substance
        assert "toxtree" in substance

        # Check regulatory structure
        assert "EU" in substance["regulatory"]
        assert "US" in substance["regulatory"]
        assert "CN" in substance["regulatory"]

        # Check toxtree structure
        assert "cramer_class" in substance["toxtree"]
        assert "cramer_value" in substance["toxtree"]
        assert "ttc" in substance["toxtree"]
        assert "cf_ttc" in substance["toxtree"]
        assert "has_alerts" in substance["toxtree"]


# ========== INTEGRATION TESTS ==========

class TestFullSessionWorkflow:
    """Test complete session workflow from creation to validation."""

    def test_complete_monolayer_workflow(self):
        """Test complete workflow for monolayer simulation."""
        # 1. Create session
        resp = client.post("/api/sessions/create", json={"name": "Monolayer Workflow"})
        session_id = resp.json()["session"]["id"]

        # 2. Add substance
        resp = client.post(
            f"/api/sessions/{session_id}/substances/unknown",
            json={"name": "additive_x", "mw": 530, "logP": 12}
        )
        substance_id = resp.json()["substance"]["id"]

        # 3. Add layer
        client.post(
            f"/api/sessions/{session_id}/layers",
            json={"polymer": "LDPE", "thickness": 100, "thickness_unit": "um"}
        )

        # 4. Add contact step (EU standard: 40°C, 10 days)
        client.post(
            f"/api/sessions/{session_id}/steps",
            json={
                "temperature": 40,
                "temperature_unit": "degC",
                "duration": 10,
                "duration_unit": "days"
            }
        )

        # 5. Assign substance to layer with C0
        client.post(
            f"/api/sessions/{session_id}/assignments",
            json={
                "substance_id": substance_id,
                "layer_index": 1,
                "C0": 1000,
                "C0_unit": "mg/kg"
            }
        )

        # 6. Validate session
        response = client.get(f"/api/sessions/{session_id}/validate")
        data = response.json()
        assert data["valid"] is True
        assert data["summary"]["layers_count"] == 1
        assert data["summary"]["substances_count"] == 1
        assert data["summary"]["assignments_count"] == 1
        assert data["summary"]["steps_count"] == 1

    def test_multilayer_setoff_workflow(self):
        """Test workflow with multilayer and set-off."""
        # 1. Create session
        resp = client.post("/api/sessions/create", json={"name": "Multilayer Set-off"})
        session_id = resp.json()["session"]["id"]

        # 2. Add substance
        resp = client.post(
            f"/api/sessions/{session_id}/substances/unknown",
            json={"name": "ink_residue", "mw": 200, "logP": 3}
        )
        substance_id = resp.json()["substance"]["id"]

        # 3. Add layers (ABA structure)
        for polymer in ["gPET", "PP", "LDPE"]:
            thickness = 20 if polymer == "gPET" else 100
            client.post(
                f"/api/sessions/{session_id}/layers",
                json={"polymer": polymer, "thickness": thickness, "thickness_unit": "um"}
            )

        # 4. Add set-off step
        client.post(
            f"/api/sessions/{session_id}/steps",
            json={
                "temperature": 25,
                "duration": 90,
                "duration_unit": "days",
                "with_food": False
            }
        )

        # 5. Add contact step
        client.post(
            f"/api/sessions/{session_id}/steps",
            json={
                "temperature": 40,
                "duration": 10,
                "duration_unit": "days",
                "with_food": True
            }
        )

        # 6. Assign substance to middle layer (PP)
        client.post(
            f"/api/sessions/{session_id}/assignments",
            json={
                "substance_id": substance_id,
                "layer_index": 2,
                "C0": 500,
                "C0_unit": "mg/kg"
            }
        )

        # 7. Validate
        response = client.get(f"/api/sessions/{session_id}/validate")
        data = response.json()
        assert data["valid"] is True
        assert data["summary"]["layers_count"] == 3
        assert data["summary"]["steps_count"] == 2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
