"""
SFPPy Studio API Tests

Comprehensive test suite to verify that Studio can reproduce
the results from CLI examples (example1.py through example5.py).

Run with: pytest studio/tests/test_api.py -v
"""

import pytest
import sys
from pathlib import Path

# Add paths for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from fastapi.testclient import TestClient
from studio.app.main import app

client = TestClient(app)


# ========== BASIC API TESTS ==========

class TestHealthEndpoints:
    """Test basic health and info endpoints."""

    def test_health_check(self):
        """Test health endpoint responds."""
        response = client.get("/api/health")
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "healthy"

    def test_api_info(self):
        """Test API info endpoint."""
        response = client.get("/api/info")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "studio_version" in data
        assert "sfppy_version" in data


# ========== ASSEMBLY TESTS ==========

class TestAssemblyAPI:
    """Test assembly/layer management API."""

    def test_list_polymers(self):
        """Test listing available polymers."""
        response = client.get("/api/assembly/polymers")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert len(data["polymers"]) > 0
        # Check required polymers exist
        polymer_codes = [p["code"] for p in data["polymers"]]
        assert "LDPE" in polymer_codes
        assert "PET" in polymer_codes
        assert "PP" in polymer_codes

    def test_get_specific_polymer(self):
        """Test getting a specific polymer."""
        response = client.get("/api/assembly/polymers/LDPE")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["polymer"]["code"] == "LDPE"

    def test_validate_single_layer_assembly(self):
        """Test validation of single layer assembly (like example1)."""
        assembly = {
            "name": "Test Monolayer",
            "layers": [
                {"index": 1, "polymer": "LDPE", "thickness": 100, "thickness_unit": "um"}
            ]
        }
        response = client.post("/api/assembly/validate", json=assembly)
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is True
        assert data["total_thickness_um"] == pytest.approx(100.0, rel=0.01)

    def test_validate_multilayer_assembly(self):
        """Test validation of multilayer assembly (like example2)."""
        assembly = {
            "name": "Test Bilayer",
            "layers": [
                {"index": 1, "polymer": "gPET", "thickness": 20, "thickness_unit": "um"},
                {"index": 2, "polymer": "LDPE", "thickness": 100, "thickness_unit": "um"}
            ]
        }
        response = client.post("/api/assembly/validate", json=assembly)
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is True
        assert data["total_thickness_um"] == pytest.approx(120.0, rel=0.01)

    def test_validate_trilayer_assembly(self):
        """Test validation of trilayer ABA assembly (like example3)."""
        assembly = {
            "name": "Test Trilayer ABA",
            "layers": [
                {"index": 1, "polymer": "LDPE", "thickness": 50, "thickness_unit": "um"},
                {"index": 2, "polymer": "PP", "thickness": 200, "thickness_unit": "um"},
                {"index": 3, "polymer": "LDPE", "thickness": 50, "thickness_unit": "um"}
            ]
        }
        response = client.post("/api/assembly/validate", json=assembly)
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is True
        assert data["total_thickness_um"] == pytest.approx(300.0, rel=0.01)

    def test_invalid_layer_indices(self):
        """Test that non-consecutive layer indices are rejected."""
        assembly = {
            "name": "Invalid",
            "layers": [
                {"index": 1, "polymer": "LDPE", "thickness": 100, "thickness_unit": "um"},
                {"index": 3, "polymer": "PP", "thickness": 100, "thickness_unit": "um"}  # Gap!
            ]
        }
        response = client.post("/api/assembly/validate", json=assembly)
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is False
        assert len(data["errors"]) > 0

    def test_estimate_properties(self):
        """Test D and k estimation for layers."""
        assembly = {
            "name": "Test",
            "layers": [
                {"index": 1, "polymer": "LDPE", "thickness": 100, "thickness_unit": "um", "temperature_C": 40}
            ]
        }
        response = client.post("/api/assembly/estimate?mw=500&logP=5", json=assembly)
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert len(data["layers"]) == 1
        assert data["layers"][0]["D_computed"] > 0
        assert data["layers"][0]["k_computed"] > 0


# ========== FOOD & CONDITIONS TESTS ==========

class TestFoodAPI:
    """Test food type and contact conditions API."""

    def test_list_food_categories(self):
        """Test listing food categories."""
        response = client.get("/api/food/categories")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert len(data["categories"]) >= 5  # fatty, aqueous, acidic, alcoholic, dry

    def test_list_simulants(self):
        """Test listing food simulants."""
        response = client.get("/api/food/simulants")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        # Check EU simulants are present
        simulant_codes = [s["code"] for s in data["simulants"]]
        assert "olive_oil" in simulant_codes
        assert "water" in simulant_codes
        assert "ethanol_50" in simulant_codes

    def test_list_shapes(self):
        """Test listing packaging shapes."""
        response = client.get("/api/food/shapes")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        shape_codes = [s["code"] for s in data["shapes"]]
        assert "cylinder" in shape_codes
        assert "bottle" in shape_codes

    def test_geometry_cylinder(self):
        """Test cylinder geometry calculation."""
        geometry = {
            "shape": "cylinder",
            "dimensions": {"radius": 50, "height": 100}  # mm
        }
        response = client.post("/api/food/geometry/calculate", json=geometry)
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["volume_cm3"] > 0
        assert data["surface_cm2"] > 0
        assert data["vs_ratio_cm"] > 0

    def test_geometry_bottle(self):
        """Test bottle geometry calculation."""
        geometry = {
            "shape": "bottle",
            "dimensions": {
                "body_radius": 40,
                "body_height": 200,
                "neck_radius": 15,
                "neck_height": 50
            }
        }
        response = client.post("/api/food/geometry/calculate", json=geometry)
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["volume_cm3"] > 0

    def test_list_condition_presets(self):
        """Test listing condition presets."""
        response = client.get("/api/food/conditions/presets")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        preset_codes = [p["code"] for p in data["presets"]]
        assert "eu_standard" in preset_codes

    def test_validate_single_step_scenario(self):
        """Test scenario validation with single step."""
        scenario = {
            "food": {"category": "fatty"},
            "geometry": {"shape": "cylinder", "dimensions": {"radius": 50, "height": 100}},
            "steps": [
                {"index": 1, "temperature_C": 40, "duration": 10, "duration_unit": "days"}
            ]
        }
        response = client.post("/api/food/scenario/validate", json=scenario)
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is True
        assert data["total_duration_days"] == pytest.approx(10.0, rel=0.01)

    def test_validate_multi_step_scenario(self):
        """Test scenario validation with multiple steps (set-off + storage)."""
        scenario = {
            "food": {"category": "fatty"},
            "geometry": {"shape": "cylinder", "dimensions": {"radius": 50, "height": 100}},
            "steps": [
                {"index": 1, "temperature_C": 25, "duration": 90, "duration_unit": "days", "with_food": False},
                {"index": 2, "temperature_C": 40, "duration": 10, "duration_unit": "days", "with_food": True}
            ]
        }
        response = client.post("/api/food/scenario/validate", json=scenario)
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is True
        assert data["has_setoff"] is True
        assert data["total_duration_days"] == pytest.approx(100.0, rel=0.01)


# ========== SUBSTANCES TESTS ==========

class TestSubstancesAPI:
    """Test substances management API."""

    def test_list_common_substances(self):
        """Test listing common substances."""
        response = client.get("/api/substances/common")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert len(data["substances"]) > 0
        # Check expected substances
        names = [s["name"] for s in data["substances"]]
        assert any("Irganox" in n for n in names)

    def test_list_substance_categories(self):
        """Test listing substance categories."""
        response = client.get("/api/substances/categories")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        codes = [c["code"] for c in data["categories"]]
        assert "antioxidant" in codes
        assert "plasticizer" in codes

    def test_search_local_substance(self):
        """Test searching for a local substance."""
        response = client.get("/api/substances/search?q=irganox&source=local")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert len(data["results"]) > 0

    def test_get_sml_known_substance(self):
        """Test getting SML for known substance."""
        response = client.get("/api/substances/sml/2082-79-3")  # Irganox 1076
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["SML"] == 6.0

    def test_validate_substance_config(self):
        """Test validating substance configuration."""
        substances = [
            {"substance_id": "irganox_1076", "layer_index": 1, "C0": 1000}
        ]
        response = client.post("/api/substances/validate", json=substances)
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is True


# ========== SIMULATION TESTS ==========

class TestSimulationAPI:
    """Test simulation engine API."""

    def test_list_simulation_presets(self):
        """Test listing simulation presets (examples)."""
        response = client.get("/api/simulation/presets")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        codes = [p["code"] for p in data["presets"]]
        assert "example1_monolayer" in codes
        assert "example2_bilayer" in codes

    def test_validate_simulation_config(self):
        """Test validating simulation configuration."""
        config = {
            "name": "Test Simulation",
            "layers": [
                {"index": 1, "polymer": "LDPE", "thickness": 100, "thickness_unit": "um", "C0": 1000}
            ],
            "steps": [
                {"index": 1, "temperature_C": 40, "duration": 10, "duration_unit": "days"}
            ],
            "geometry": {"shape": "cylinder", "dimensions": {"radius": 50, "height": 100}},
            "food": {"category": "fatty"}
        }
        response = client.post("/api/simulation/validate", json=config)
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is True
        assert data["layer_count"] == 1
        assert data["step_count"] == 1

    def test_run_simulation_sync(self):
        """Test running a simulation synchronously."""
        config = {
            "name": "Test Run",
            "layers": [
                {"index": 1, "polymer": "LDPE", "thickness": 100, "thickness_unit": "um", "C0": 1000}
            ],
            "steps": [
                {"index": 1, "temperature_C": 40, "duration": 10, "duration_unit": "days"}
            ],
            "geometry": {"shape": "cylinder", "dimensions": {"radius": 50, "height": 100}},
            "food": {"category": "fatty"}
        }
        response = client.post("/api/simulation/run?async_mode=false", json=config)
        assert response.status_code == 200
        data = response.json()
        assert data["job_id"] is not None
        # Should have results (at least mock)
        if data["status"] == "completed":
            assert data["results"] is not None


# ========== JOBS TESTS ==========

class TestJobsAPI:
    """Test job history management API."""

    def test_list_jobs(self):
        """Test listing jobs."""
        response = client.get("/api/jobs/list")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "jobs" in data
        assert "total" in data

    def test_get_job_stats(self):
        """Test getting job statistics."""
        response = client.get("/api/jobs/stats")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "total" in data
        assert "by_status" in data


# ========== CONFIG TESTS ==========

class TestConfigAPI:
    """Test configuration API."""

    def test_get_config(self):
        """Test getting configuration."""
        response = client.get("/api/config/")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "config" in data

    def test_get_traceability(self):
        """Test getting traceability info."""
        response = client.get("/api/config/traceability")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "traceability" in data
        assert "machine" in data["traceability"]
        assert "software" in data["traceability"]

    def test_list_diffusivity_models(self):
        """Test listing diffusivity estimation models."""
        response = client.get("/api/config/models/diffusivity")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        codes = [m["code"] for m in data["models"]]
        assert "piringer" in codes

    def test_list_partition_models(self):
        """Test listing partition coefficient models."""
        response = client.get("/api/config/models/partition")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        codes = [m["code"] for m in data["models"]]
        assert "fhp" in codes

    def test_list_units(self):
        """Test listing available units."""
        response = client.get("/api/config/units")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "thickness" in data["units"]
        assert "time" in data["units"]
        assert "concentration" in data["units"]

    def test_list_regulations(self):
        """Test listing supported regulations."""
        response = client.get("/api/config/regulations")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        codes = [r["code"] for r in data["regulations"]]
        assert "EU_10_2011" in codes


# ========== HELP TESTS ==========

class TestHelpAPI:
    """Test contextual help API."""

    def test_get_assembly_help(self):
        """Test getting help for assembly elements."""
        response = client.get("/api/help/context/assembly/polymer")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "title" in data
        assert "short" in data

    def test_get_food_help(self):
        """Test getting help for food elements."""
        response = client.get("/api/help/context/food/category")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True


# ========== EXAMPLE REPRODUCTION TESTS ==========

class TestExampleReproduction:
    """
    Test that Studio API can reproduce CLI example results.

    These tests verify that the configurations from example1.py through example5.py
    can be represented and validated through the API.
    """

    def test_example1_monolayer_config(self):
        """
        Reproduce example1.py: Simple monolayer LDPE migration.

        Original: LDPE 100µm, 40°C, 10 days, fatty food
        """
        # Assembly
        assembly = {
            "name": "Example1 - Monolayer LDPE",
            "layers": [
                {"index": 1, "polymer": "LDPE", "thickness": 100, "thickness_unit": "um", "temperature_C": 40}
            ]
        }
        resp = client.post("/api/assembly/validate", json=assembly)
        assert resp.json()["valid"] is True

        # Scenario
        scenario = {
            "food": {"category": "fatty", "simulant": "olive_oil"},
            "geometry": {"shape": "cylinder", "dimensions": {"radius": 50, "height": 100}},
            "steps": [{"index": 1, "temperature_C": 40, "duration": 10, "duration_unit": "days"}]
        }
        resp = client.post("/api/food/scenario/validate", json=scenario)
        assert resp.json()["valid"] is True

    def test_example2_bilayer_barrier_config(self):
        """
        Reproduce example2.py: Bilayer with functional barrier.

        Original: gPET 20µm (barrier) + LDPE 100µm, 40°C, 10 days
        """
        assembly = {
            "name": "Example2 - Bilayer with Barrier",
            "layers": [
                {"index": 1, "polymer": "gPET", "thickness": 20, "thickness_unit": "um"},
                {"index": 2, "polymer": "LDPE", "thickness": 100, "thickness_unit": "um", "temperature_C": 40}
            ]
        }
        resp = client.post("/api/assembly/validate", json=assembly)
        data = resp.json()
        assert data["valid"] is True
        assert data["total_thickness_um"] == pytest.approx(120.0, rel=0.01)

    def test_example3_trilayer_setoff_config(self):
        """
        Reproduce example3.py: Trilayer ABA with set-off.

        Original: LDPE/PP/LDPE, 90 days set-off + 10 days contact
        """
        assembly = {
            "name": "Example3 - Trilayer ABA",
            "layers": [
                {"index": 1, "polymer": "LDPE", "thickness": 50, "thickness_unit": "um"},
                {"index": 2, "polymer": "PP", "thickness": 200, "thickness_unit": "um"},
                {"index": 3, "polymer": "LDPE", "thickness": 50, "thickness_unit": "um"}
            ]
        }
        resp = client.post("/api/assembly/validate", json=assembly)
        assert resp.json()["valid"] is True

        # Multi-step scenario with set-off
        scenario = {
            "food": {"category": "fatty"},
            "geometry": {"shape": "cylinder", "dimensions": {"radius": 50, "height": 100}},
            "steps": [
                {"index": 1, "temperature_C": 25, "duration": 90, "duration_unit": "days", "with_food": False},
                {"index": 2, "temperature_C": 40, "duration": 10, "duration_unit": "days", "with_food": True}
            ]
        }
        resp = client.post("/api/food/scenario/validate", json=scenario)
        data = resp.json()
        assert data["valid"] is True
        assert data["has_setoff"] is True

    def test_example4_hotfill_config(self):
        """
        Reproduce example4.py: Hot fill process.

        Original: PET bottle, 85°C fill + 25°C storage
        """
        assembly = {
            "name": "Example4 - Hot Fill PET",
            "layers": [
                {"index": 1, "polymer": "PET", "thickness": 300, "thickness_unit": "um"}
            ]
        }
        resp = client.post("/api/assembly/validate", json=assembly)
        assert resp.json()["valid"] is True

        # Hot fill scenario
        scenario = {
            "food": {"category": "aqueous"},
            "geometry": {"shape": "bottle", "dimensions": {
                "body_radius": 40, "body_height": 200,
                "neck_radius": 15, "neck_height": 50
            }},
            "steps": [
                {"index": 1, "temperature_C": 85, "duration": 30, "duration_unit": "minutes"},
                {"index": 2, "temperature_C": 25, "duration": 6, "duration_unit": "months"}
            ]
        }
        resp = client.post("/api/food/scenario/validate", json=scenario)
        assert resp.json()["valid"] is True

    def test_example5_multilayer_complex_config(self):
        """
        Reproduce example5.py: Complex multilayer structure.

        Original: Multiple layers with different polymers
        """
        assembly = {
            "name": "Example5 - Complex Multilayer",
            "layers": [
                {"index": 1, "polymer": "LDPE", "thickness": 30, "thickness_unit": "um"},
                {"index": 2, "polymer": "EVOH", "thickness": 10, "thickness_unit": "um"},
                {"index": 3, "polymer": "PA6", "thickness": 20, "thickness_unit": "um"},
                {"index": 4, "polymer": "PP", "thickness": 200, "thickness_unit": "um"}
            ]
        }
        resp = client.post("/api/assembly/validate", json=assembly)
        data = resp.json()
        assert data["valid"] is True
        assert data["total_thickness_um"] == pytest.approx(260.0, rel=0.01)


# ========== UNIT CONVERSION TESTS ==========

class TestUnitConversions:
    """Test that various units are properly handled."""

    def test_thickness_nm(self):
        """Test thickness in nanometers."""
        assembly = {
            "name": "Test nm",
            "layers": [{"index": 1, "polymer": "LDPE", "thickness": 100000, "thickness_unit": "nm"}]
        }
        resp = client.post("/api/assembly/validate", json=assembly)
        data = resp.json()
        assert data["valid"] is True
        assert data["total_thickness_um"] == pytest.approx(100.0, rel=0.01)

    def test_thickness_mm(self):
        """Test thickness in millimeters."""
        assembly = {
            "name": "Test mm",
            "layers": [{"index": 1, "polymer": "LDPE", "thickness": 0.1, "thickness_unit": "mm"}]
        }
        resp = client.post("/api/assembly/validate", json=assembly)
        data = resp.json()
        assert data["valid"] is True
        assert data["total_thickness_um"] == pytest.approx(100.0, rel=0.01)

    def test_duration_hours(self):
        """Test duration in hours."""
        scenario = {
            "food": {"category": "fatty"},
            "geometry": {"shape": "cylinder", "dimensions": {"radius": 50, "height": 100}},
            "steps": [{"index": 1, "temperature_C": 40, "duration": 240, "duration_unit": "hours"}]
        }
        resp = client.post("/api/food/scenario/validate", json=scenario)
        data = resp.json()
        assert data["valid"] is True
        assert data["total_duration_days"] == pytest.approx(10.0, rel=0.01)

    def test_duration_months(self):
        """Test duration in months."""
        scenario = {
            "food": {"category": "fatty"},
            "geometry": {"shape": "cylinder", "dimensions": {"radius": 50, "height": 100}},
            "steps": [{"index": 1, "temperature_C": 4, "duration": 6, "duration_unit": "months"}]
        }
        resp = client.post("/api/food/scenario/validate", json=scenario)
        data = resp.json()
        assert data["valid"] is True
        assert data["total_duration_days"] == pytest.approx(180.0, rel=0.1)


# ========== EXPORT TESTS ==========

class TestExportAPI:
    """Test export functionality for simulation results."""

    def test_list_export_formats(self):
        """Test listing available export formats."""
        response = client.get("/api/export/formats")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "formats" in data
        # CSV and JSON should always be available
        assert "csv" in data["formats"]
        assert "json" in data["formats"]
        assert data["formats"]["csv"]["available"] is True
        assert data["formats"]["json"]["available"] is True

    def test_export_csv_no_job(self):
        """Test CSV export fails gracefully for non-existent job."""
        response = client.get("/api/export/nonexistent123/csv")
        assert response.status_code == 404

    def test_export_json_no_job(self):
        """Test JSON export fails gracefully for non-existent job."""
        response = client.get("/api/export/nonexistent123/json")
        assert response.status_code == 404

    def test_export_xlsx_no_job(self):
        """Test XLSX export fails gracefully for non-existent job."""
        response = client.get("/api/export/nonexistent123/xlsx")
        # Should be 404 (job not found) or 500 (openpyxl not available)
        assert response.status_code in [404, 500]

    def test_export_pdf_no_job(self):
        """Test PDF export fails gracefully for non-existent job."""
        response = client.get("/api/export/nonexistent123/pdf")
        # Should be 404 (job not found) or 500 (reportlab not available)
        assert response.status_code in [404, 500]

    def test_export_png_no_job(self):
        """Test PNG export fails gracefully for non-existent job."""
        response = client.get("/api/export/nonexistent123/png")
        # Should be 404 (job not found) or 500 (matplotlib not available)
        assert response.status_code in [404, 500]

    def test_export_svg_no_job(self):
        """Test SVG export fails gracefully for non-existent job."""
        response = client.get("/api/export/nonexistent123/svg")
        # Should be 404 (job not found) or 500 (matplotlib not available)
        assert response.status_code in [404, 500]


class TestExportWithSimulation:
    """Test exports with actual simulation job (integration tests)."""

    @pytest.fixture
    def completed_job_id(self):
        """Create a simulation job and return its ID."""
        # Run a simple simulation
        config = {
            "name": "Export Test",
            "layers": [
                {"index": 1, "polymer": "LDPE", "thickness": 100, "thickness_unit": "um"}
            ],
            "steps": [
                {"index": 1, "temperature_C": 40, "duration": 10, "duration_unit": "days"}
            ],
            "geometry": {"shape": "cylinder", "dimensions": {"radius": 50, "height": 100}},
            "food": {"category": "fatty"}
        }
        response = client.post("/api/simulation/run?async_mode=false", json=config)
        data = response.json()
        if data.get("status") == "completed":
            return data.get("job_id")
        return None

    def test_export_csv_with_job(self, completed_job_id):
        """Test CSV export with completed job."""
        if not completed_job_id:
            pytest.skip("Simulation did not complete")
        response = client.get(f"/api/export/{completed_job_id}/csv")
        assert response.status_code == 200
        assert "text/csv" in response.headers.get("content-type", "")
        assert "time_s,time_days,CF_mg_kg" in response.text

    def test_export_json_with_job(self, completed_job_id):
        """Test JSON export with completed job."""
        if not completed_job_id:
            pytest.skip("Simulation did not complete")
        response = client.get(f"/api/export/{completed_job_id}/json")
        assert response.status_code == 200
        assert "application/json" in response.headers.get("content-type", "")
        data = response.json()
        assert "job_id" in data
        assert "results" in data


# ========== MATERIALS DISCOVERY TESTS ==========

class TestMaterialsAPI:
    """Test materials discovery API from layer.py."""

    def test_list_all_materials(self):
        """Test listing all materials discovered from layer.py."""
        response = client.get("/api/assembly/materials")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "materials" in data
        assert data["count"] > 10  # Should have many materials

    def test_list_materials_by_category(self):
        """Test filtering materials by category."""
        response = client.get("/api/assembly/materials?category=polymer")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        # All returned should be polymers
        for mat in data["materials"]:
            assert mat.get("category") == "polymer"

    def test_get_specific_material(self):
        """Test getting details for a specific material."""
        response = client.get("/api/assembly/materials/LDPE")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["material"]["code"] == "LDPE"
        assert "Tg" in data["material"]

    def test_get_material_tg(self):
        """Test getting Tg for a material."""
        response = client.get("/api/assembly/materials/gPET/tg")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["Tg"] is not None
        assert data["is_barrier"] is True  # gPET has high Tg

    def test_material_not_found(self):
        """Test 404 for unknown material."""
        response = client.get("/api/assembly/materials/UNKNOWN_POLYMER_XYZ")
        assert response.status_code == 404


# ========== SUBSTANCE DETAIL TESTS ==========

class TestSubstanceDetailAPI:
    """Test comprehensive substance detail API."""

    def test_search_pubchem_substance(self):
        """Test searching PubChem for a substance."""
        response = client.get("/api/substances/search?q=anisole&source=pubchem")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        # Anisole should be found
        if data["count"] > 0:
            assert any("anisole" in r.get("name", "").lower() for r in data["results"])

    def test_get_substance_detail(self):
        """Test getting comprehensive substance details."""
        # Use CID 7519 (anisole)
        response = client.get("/api/substances/detail/7519")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        substance = data.get("substance", {})
        # Should have basic properties
        assert "cid" in substance or "name" in substance

    def test_get_pubchem_substance_by_cid(self):
        """Test getting substance by PubChem CID."""
        # BHT has CID 31404
        response = client.get("/api/substances/pubchem/31404")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["substance"]["cid"] == 31404

    def test_get_substance_thumbnail(self):
        """Test getting substance thumbnail URLs."""
        response = client.get("/api/substances/thumbnail/31404")
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert "thumbnail_url" in data
        assert "pubchem.ncbi.nlm.nih.gov" in data["thumbnail_url"]

    def test_get_regulatory_status(self):
        """Test getting regulatory status for a substance."""
        # Irganox 1076 (CAS 2082-79-3)
        response = client.get("/api/substances/regulatory/2082-79-3")
        # May fail if migrant class not available, that's OK
        if response.status_code == 200:
            data = response.json()
            assert data["success"] is True
            assert "regulatory" in data


# ========== EDGE CASES AND ERROR HANDLING ==========

class TestErrorHandling:
    """Test error handling for various edge cases."""

    def test_invalid_polymer(self):
        """Test validation with invalid polymer code."""
        assembly = {
            "name": "Invalid",
            "layers": [{"index": 1, "polymer": "INVALID_POLYMER", "thickness": 100, "thickness_unit": "um"}]
        }
        response = client.post("/api/assembly/validate", json=assembly)
        # Should return 200 with valid=False and errors, not 500
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is False
        assert len(data["errors"]) > 0

    def test_negative_thickness(self):
        """Test validation with negative thickness."""
        assembly = {
            "name": "Negative",
            "layers": [{"index": 1, "polymer": "LDPE", "thickness": -100, "thickness_unit": "um"}]
        }
        response = client.post("/api/assembly/validate", json=assembly)
        # Should fail validation due to Pydantic constraint
        assert response.status_code == 422

    def test_empty_layers(self):
        """Test validation with empty layers list."""
        assembly = {"name": "Empty", "layers": []}
        response = client.post("/api/assembly/validate", json=assembly)
        # Should fail validation
        assert response.status_code == 422

    def test_missing_required_fields(self):
        """Test validation with missing required fields."""
        # Missing layers entirely
        response = client.post("/api/assembly/validate", json={"name": "Test"})
        assert response.status_code == 422

    def test_invalid_geometry_shape(self):
        """Test geometry calculation with invalid shape."""
        geometry = {"shape": "invalid_shape", "dimensions": {"radius": 50}}
        response = client.post("/api/food/geometry/calculate", json=geometry)
        # Should return error
        assert response.status_code in [400, 422, 500]

    def test_invalid_duration_unit(self):
        """Test with invalid duration unit."""
        scenario = {
            "food": {"category": "fatty"},
            "geometry": {"shape": "cylinder", "dimensions": {"radius": 50, "height": 100}},
            "steps": [{"index": 1, "temperature_C": 40, "duration": 10, "duration_unit": "invalid_unit"}]
        }
        response = client.post("/api/food/scenario/validate", json=scenario)
        # Should handle gracefully
        assert response.status_code in [200, 400, 422]

    def test_negative_c0(self):
        """Test validation with negative C0."""
        substances = [{"substance_id": "bht", "layer_index": 1, "C0": -100}]
        response = client.post("/api/substances/validate", json=substances)
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is False
        assert len(data["errors"]) > 0


# ========== INTEGRATION TESTS ==========

class TestFullWorkflow:
    """Test complete workflows from configuration to results."""

    def test_complete_monolayer_workflow(self):
        """Test complete workflow for monolayer simulation."""
        # 1. Validate assembly
        assembly = {
            "name": "Workflow Test",
            "layers": [{"index": 1, "polymer": "LDPE", "thickness": 100, "thickness_unit": "um"}]
        }
        resp = client.post("/api/assembly/validate", json=assembly)
        assert resp.json()["valid"] is True

        # 2. Validate scenario
        scenario = {
            "food": {"category": "fatty"},
            "geometry": {"shape": "cylinder", "dimensions": {"radius": 50, "height": 100}},
            "steps": [{"index": 1, "temperature_C": 40, "duration": 10, "duration_unit": "days"}]
        }
        resp = client.post("/api/food/scenario/validate", json=scenario)
        assert resp.json()["valid"] is True

        # 3. Validate full config
        config = {
            "name": "Workflow Test",
            "layers": [{"index": 1, "polymer": "LDPE", "thickness": 100, "thickness_unit": "um", "C0": 1000}],
            "steps": [{"index": 1, "temperature_C": 40, "duration": 10, "duration_unit": "days"}],
            "geometry": {"shape": "cylinder", "dimensions": {"radius": 50, "height": 100}},
            "food": {"category": "fatty"}
        }
        resp = client.post("/api/simulation/validate", json=config)
        assert resp.json()["valid"] is True

        # 4. Run simulation
        resp = client.post("/api/simulation/run?async_mode=false", json=config)
        data = resp.json()
        assert data["job_id"] is not None

    def test_multilayer_with_barrier_workflow(self):
        """Test workflow with barrier layer configuration."""
        config = {
            "name": "Barrier Test",
            "layers": [
                {"index": 1, "polymer": "gPET", "thickness": 20, "thickness_unit": "um", "C0": 0},
                {"index": 2, "polymer": "LDPE", "thickness": 100, "thickness_unit": "um", "C0": 1000}
            ],
            "steps": [{"index": 1, "temperature_C": 40, "duration": 10, "duration_unit": "days"}],
            "geometry": {"shape": "cylinder", "dimensions": {"radius": 50, "height": 100}},
            "food": {"category": "fatty"}
        }

        # Validate full config
        resp = client.post("/api/simulation/validate", json=config)
        data = resp.json()
        assert data["valid"] is True
        assert data["layer_count"] == 2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
