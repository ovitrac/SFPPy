"""
Tests for SFPPy Studio Session Models

Phase 1 tests: Validate Pydantic models and example JSON files.
Ensures support for multiple substances, layers, and contact steps.

@author: SFPPy Studio
@license: MIT
"""

import pytest
import json
from pathlib import Path
from datetime import datetime

# Import session models
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.models.session import (
    Session,
    Metadata,
    Geometry,
    GeometryShape,
    GeometryDimensions,
    Substance,
    SubstanceSource,
    SubstanceProperties,
    Layer,
    LayerSubstance,
    PolymerType,
    Food,
    FoodCategory,
    FoodTexture,
    FoodAffinity,
    ContactStep,
    ContactType,
    SimulationConfig,
    ValueWithUnit,
    create_empty_session,
)


# =============================================================================
# Test Fixtures
# =============================================================================

@pytest.fixture
def examples_dir():
    """Path to example JSON files."""
    return Path(__file__).parent.parent / "examples"


@pytest.fixture
def example1_path(examples_dir):
    """Path to example1.sfppy.json."""
    return examples_dir / "example1.sfppy.json"


@pytest.fixture
def example1_extensions_path(examples_dir):
    """Path to example1_extensions.sfppy.json."""
    return examples_dir / "example1_extensions.sfppy.json"


@pytest.fixture
def example3_path(examples_dir):
    """Path to example3.sfppy.json."""
    return examples_dir / "example3.sfppy.json"


@pytest.fixture
def example3_variant_path(examples_dir):
    """Path to example3_variant.sfppy.json."""
    return examples_dir / "example3_variant.sfppy.json"


# =============================================================================
# Phase 1 Tests: Basic Model Validation
# =============================================================================

class TestValueWithUnit:
    """Test ValueWithUnit model."""

    def test_create_value_with_unit(self):
        """Test creating a value with unit."""
        v = ValueWithUnit(value=100, unit="um")
        assert v.value == 100
        assert v.unit == "um"

    def test_value_serialization(self):
        """Test JSON serialization."""
        v = ValueWithUnit(value=25.5, unit="degC")
        data = v.model_dump()
        assert data == {"value": 25.5, "unit": "degC"}


class TestGeometry:
    """Test Geometry models."""

    def test_cylinder_geometry(self):
        """Test cylinder geometry creation."""
        geom = Geometry(
            shape=GeometryShape.CYLINDER,
            dimensions=GeometryDimensions(
                radius=ValueWithUnit(value=30, unit="mm"),
                height=ValueWithUnit(value=19, unit="cm")
            )
        )
        assert geom.shape == GeometryShape.CYLINDER
        assert geom.dimensions.radius.value == 30
        assert geom.dimensions.height.value == 19

    def test_box_container_geometry(self):
        """Test box container geometry creation."""
        geom = Geometry(
            shape=GeometryShape.BOX_CONTAINER,
            dimensions=GeometryDimensions(
                length=ValueWithUnit(value=19, unit="cm"),
                width=ValueWithUnit(value=10, unit="cm"),
                height=ValueWithUnit(value=8, unit="cm")
            )
        )
        assert geom.shape == GeometryShape.BOX_CONTAINER
        assert geom.dimensions.length.value == 19


class TestSubstance:
    """Test Substance models."""

    def test_create_pubchem_substance(self):
        """Test creating a PubChem substance."""
        substance = Substance(
            id="m1",
            source=SubstanceSource.PUBCHEM,
            lookup_name="irganox 1076",
            cid=16386,
            properties=SubstanceProperties(
                name="Irganox 1076",
                cas=["2082-79-3"],
                mw=530.9,
                logP=13.8
            ),
            SML=6.0
        )
        assert substance.id == "m1"
        assert substance.cid == 16386
        assert substance.properties.mw == 530.9
        assert substance.SML == 6.0

    def test_substance_auto_id(self):
        """Test that substance gets auto-generated ID."""
        substance = Substance(
            source=SubstanceSource.MANUAL,
            properties=SubstanceProperties(name="Test")
        )
        assert substance.id is not None
        assert len(substance.id) == 8


class TestLayer:
    """Test Layer models."""

    def test_create_layer(self):
        """Test creating a polymer layer."""
        layer = Layer(
            index=1,
            polymer=PolymerType.LDPE,
            thickness=ValueWithUnit(value=100, unit="um"),
            name="LDPE film"
        )
        assert layer.index == 1
        assert layer.polymer == PolymerType.LDPE
        assert layer.thickness.value == 100

    def test_layer_with_substance(self):
        """Test layer with substance assignment."""
        layer = Layer(
            index=1,
            polymer=PolymerType.PP,
            thickness=ValueWithUnit(value=500, unit="um"),
            substances=[
                LayerSubstance(
                    substance_id="limonene",
                    C0=ValueWithUnit(value=200, unit="mg/kg"),
                    D_auto=True,
                    k_auto=True
                )
            ]
        )
        assert len(layer.substances) == 1
        assert layer.substances[0].C0.value == 200

    def test_layer_index_validation(self):
        """Test that layer index must be >= 1."""
        with pytest.raises(ValueError):
            Layer(
                index=0,  # Invalid: must be >= 1
                polymer=PolymerType.LDPE,
                thickness=ValueWithUnit(value=100, unit="um")
            )


class TestContactStep:
    """Test ContactStep models."""

    def test_create_storage_step(self):
        """Test creating a storage contact step."""
        step = ContactStep(
            index=1,
            name="Cold storage",
            type=ContactType.STORAGE,
            temperature=ValueWithUnit(value=7, unit="degC"),
            duration=ValueWithUnit(value=10, unit="days")
        )
        assert step.type == ContactType.STORAGE
        assert step.temperature.value == 7
        assert step.with_food_contact == True  # Default

    def test_setoff_step(self):
        """Test creating a setoff step (no food contact)."""
        step = ContactStep(
            index=1,
            name="Setoff storage",
            type=ContactType.SETOFF,
            temperature=ValueWithUnit(value=20, unit="degC"),
            duration=ValueWithUnit(value=4, unit="months"),
            with_food_contact=False
        )
        assert step.type == ContactType.SETOFF
        assert step.with_food_contact == False


class TestFood:
    """Test Food models."""

    def test_create_realfood(self):
        """Test creating a real food definition."""
        food = Food(
            category=FoodCategory.REALFOOD,
            texture=FoodTexture.SEMISOLID,
            affinity=FoodAffinity.FAT,
            simulant="ethanol",
            name="Fatty sandwich"
        )
        assert food.category == FoodCategory.REALFOOD
        assert food.affinity == FoodAffinity.FAT


# =============================================================================
# Phase 1 Tests: Session Validation
# =============================================================================

class TestSession:
    """Test complete Session models."""

    def test_create_empty_session(self):
        """Test creating an empty session with defaults."""
        session = create_empty_session("Test Session")
        assert session.metadata.name == "Test Session"
        assert len(session.layers) == 1
        assert len(session.contact_steps) == 1

    def test_session_requires_layers(self):
        """Test that session requires at least one layer."""
        with pytest.raises(ValueError):
            Session(
                metadata=Metadata(name="Invalid"),
                geometry=Geometry(
                    shape=GeometryShape.CYLINDER,
                    dimensions=GeometryDimensions(
                        radius=ValueWithUnit(value=30, unit="mm"),
                        height=ValueWithUnit(value=19, unit="cm")
                    )
                ),
                layers=[],  # Empty - should fail
                contact_steps=[
                    ContactStep(
                        index=1,
                        temperature=ValueWithUnit(value=25, unit="degC"),
                        duration=ValueWithUnit(value=1, unit="days")
                    )
                ]
            )

    def test_session_requires_contact_steps(self):
        """Test that session requires at least one contact step."""
        with pytest.raises(ValueError):
            Session(
                metadata=Metadata(name="Invalid"),
                geometry=Geometry(
                    shape=GeometryShape.CYLINDER,
                    dimensions=GeometryDimensions(
                        radius=ValueWithUnit(value=30, unit="mm"),
                        height=ValueWithUnit(value=19, unit="cm")
                    )
                ),
                layers=[
                    Layer(
                        index=1,
                        polymer=PolymerType.LDPE,
                        thickness=ValueWithUnit(value=100, unit="um")
                    )
                ],
                contact_steps=[]  # Empty - should fail
            )


# =============================================================================
# Phase 1 Tests: Example JSON Validation
# =============================================================================

class TestExampleFiles:
    """Test loading and validating example JSON files."""

    def test_example1_loads(self, example1_path):
        """Test that example1.sfppy.json loads correctly."""
        assert example1_path.exists(), f"Example file not found: {example1_path}"
        session = Session.load(str(example1_path))
        assert session.metadata.name == "Example 1: LDPE Film Migration"

    def test_example1_multiple_substances(self, example1_path):
        """Test that example1 has multiple substances."""
        session = Session.load(str(example1_path))
        assert len(session.substances) == 2
        substance_ids = [s.id for s in session.substances]
        assert "m1" in substance_ids
        assert "m2" in substance_ids

    def test_example1_substance_properties(self, example1_path):
        """Test substance properties in example1."""
        session = Session.load(str(example1_path))
        m1 = next(s for s in session.substances if s.id == "m1")
        assert m1.properties.name == "Irganox 1076"
        assert m1.properties.mw == 530.9
        assert m1.SML == 6.0

    def test_example1_extensions_multiple_steps(self, example1_extensions_path):
        """Test that example1_extensions has multiple contact steps."""
        session = Session.load(str(example1_extensions_path))
        assert len(session.contact_steps) == 2
        assert session.contact_steps[0].name == "Cold storage"
        assert session.contact_steps[1].name == "Warming before consumption"

    def test_example3_loads(self, example3_path):
        """Test that example3.sfppy.json loads correctly."""
        session = Session.load(str(example3_path))
        assert session.metadata.name == "Example 3: Trilayer ABA Migration"

    def test_example3_multiple_layers(self, example3_path):
        """Test that example3 has multiple layers (ABA structure)."""
        session = Session.load(str(example3_path))
        assert len(session.layers) == 3
        # Check ABA structure: wPET / PP / gPET
        assert session.layers[0].polymer == PolymerType.wPET
        assert session.layers[1].polymer == PolymerType.PP
        assert session.layers[2].polymer == PolymerType.gPET

    def test_example3_layer_thicknesses(self, example3_path):
        """Test layer thicknesses in example3."""
        session = Session.load(str(example3_path))
        # wPET: 20 um, PP: 500 um, gPET: 20 um
        assert session.layers[0].thickness.value == 20
        assert session.layers[1].thickness.value == 500
        assert session.layers[2].thickness.value == 20

    def test_example3_multiple_contact_steps(self, example3_path):
        """Test that example3 has three contact steps."""
        session = Session.load(str(example3_path))
        assert len(session.contact_steps) == 3
        # Check step types
        assert session.contact_steps[0].type == ContactType.SETOFF
        assert session.contact_steps[1].type == ContactType.HOTFILLING
        assert session.contact_steps[2].type == ContactType.STORAGE

    def test_example3_setoff_no_food_contact(self, example3_path):
        """Test that setoff step has no food contact."""
        session = Session.load(str(example3_path))
        setoff_step = session.contact_steps[0]
        assert setoff_step.with_food_contact == False

    def test_example3_variant_reduced_thickness(self, example3_variant_path):
        """Test variant with reduced barrier thickness."""
        session = Session.load(str(example3_variant_path))
        # PET layers should be 10 um instead of 20 um
        assert session.layers[0].thickness.value == 10
        assert session.layers[2].thickness.value == 10


# =============================================================================
# Phase 1 Tests: Serialization Roundtrip
# =============================================================================

class TestSerialization:
    """Test JSON serialization and deserialization."""

    def test_roundtrip_example1(self, example1_path):
        """Test that example1 survives JSON roundtrip."""
        session = Session.load(str(example1_path))
        json_str = session.to_json()
        reloaded = Session.from_json(json_str)

        assert reloaded.metadata.name == session.metadata.name
        assert len(reloaded.substances) == len(session.substances)
        assert len(reloaded.layers) == len(session.layers)

    def test_roundtrip_example3(self, example3_path):
        """Test that example3 survives JSON roundtrip."""
        session = Session.load(str(example3_path))
        json_str = session.to_json()
        reloaded = Session.from_json(json_str)

        assert reloaded.metadata.name == session.metadata.name
        assert len(reloaded.layers) == 3
        assert len(reloaded.contact_steps) == 3

    def test_roundtrip_preserves_substance_assignments(self, example3_path):
        """Test that substance assignments are preserved."""
        session = Session.load(str(example3_path))
        json_str = session.to_json()
        reloaded = Session.from_json(json_str)

        # Check PP layer has limonene at 200 mg/kg
        pp_layer = reloaded.layers[1]
        assert len(pp_layer.substances) == 1
        assert pp_layer.substances[0].C0.value == 200


# =============================================================================
# Phase 1 Tests: Edge Cases
# =============================================================================

class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_layer_multiple_substances(self):
        """Test layer with multiple substances."""
        layer = Layer(
            index=1,
            polymer=PolymerType.LDPE,
            thickness=ValueWithUnit(value=100, unit="um"),
            substances=[
                LayerSubstance(
                    substance_id="m1",
                    C0=ValueWithUnit(value=5000, unit="mg/kg")
                ),
                LayerSubstance(
                    substance_id="m2",
                    C0=ValueWithUnit(value=5000, unit="mg/kg")
                )
            ]
        )
        assert len(layer.substances) == 2

    def test_substance_with_override_d(self):
        """Test substance with manual D override."""
        layer_substance = LayerSubstance(
            substance_id="test",
            C0=ValueWithUnit(value=1000, unit="mg/kg"),
            D_auto=False,
            D_override=1.5e-13
        )
        assert layer_substance.D_auto == False
        assert layer_substance.D_override == 1.5e-13

    def test_all_polymer_types(self):
        """Test that all polymer types are valid."""
        for ptype in PolymerType:
            layer = Layer(
                index=1,
                polymer=ptype,
                thickness=ValueWithUnit(value=100, unit="um")
            )
            assert layer.polymer == ptype

    def test_all_geometry_shapes(self):
        """Test that all geometry shapes are valid."""
        for shape in GeometryShape:
            geom = Geometry(
                shape=shape,
                dimensions=GeometryDimensions(
                    radius=ValueWithUnit(value=30, unit="mm"),
                    height=ValueWithUnit(value=10, unit="cm")
                )
            )
            assert geom.shape == shape


# =============================================================================
# Run tests
# =============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
