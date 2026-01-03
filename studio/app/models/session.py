"""
Session Data Models for SFPPy Studio

Defines Pydantic models for the .sfppy.json session format.
Supports save/load of complete simulation configurations.

Version: 1.0
"""

from datetime import datetime
from typing import Optional, List, Dict, Any, Union
from pydantic import BaseModel, Field, field_validator
from enum import Enum
import uuid


# =============================================================================
# Enums
# =============================================================================

class GeometryShape(str, Enum):
    CYLINDER = "cylinder"
    BOTTLE = "bottle"
    BOX_CONTAINER = "box_container"
    POUCH = "pouch"
    SPHERE = "sphere"


class PolymerType(str, Enum):
    LDPE = "LDPE"
    LLDPE = "LLDPE"
    HDPE = "HDPE"
    PP = "PP"
    gPET = "gPET"
    wPET = "wPET"
    rPET = "rPET"
    PS = "PS"
    PMMA = "PMMA"
    PA6 = "PA6"
    PA66 = "PA66"
    PBT = "PBT"
    PEN = "PEN"
    HIPS = "HIPS"
    SBS = "SBS"


class FoodCategory(str, Enum):
    REALFOOD = "realfood"
    SIMULANT = "simulant"
    SETOFF = "setoff"


class FoodTexture(str, Enum):
    LIQUID = "liquid"
    SEMISOLID = "semisolid"
    SOLID = "solid"


class FoodAffinity(str, Enum):
    FAT = "fat"
    AQUEOUS = "aqueous"
    ACIDIC = "acidic"
    ALCOHOLIC = "alcoholic"


class ContactType(str, Enum):
    STORAGE = "storage"
    HOTFILLING = "hotfilling"
    SETOFF = "setoff"
    AMBIENT = "ambient"


class SubstanceSource(str, Enum):
    PUBCHEM = "pubchem"
    MANUAL = "manual"
    INTERNAL = "internal"


# =============================================================================
# Value with Unit Models
# =============================================================================

class ValueWithUnit(BaseModel):
    """A numeric value with explicit unit."""
    value: float
    unit: str

    def to_si(self) -> float:
        """Convert to SI units."""
        from studio.app.utils.units import convert_to_si
        return convert_to_si(self.value, self.unit)


class OptionalValueWithUnit(BaseModel):
    """Optional value with unit - can be auto-computed."""
    value: Optional[float] = None
    unit: Optional[str] = None
    auto: bool = True


# =============================================================================
# Geometry Models
# =============================================================================

class GeometryDimensions(BaseModel):
    """Geometry dimensions - shape-specific fields."""
    # Cylinder
    radius: Optional[ValueWithUnit] = None
    height: Optional[ValueWithUnit] = None
    length: Optional[ValueWithUnit] = None  # Alternative to height for cylinder

    # Bottle
    body_radius: Optional[ValueWithUnit] = None
    body_height: Optional[ValueWithUnit] = None
    neck_radius: Optional[ValueWithUnit] = None
    neck_height: Optional[ValueWithUnit] = None

    # Box/Rectangle
    width: Optional[ValueWithUnit] = None

    # Pouch
    thickness: Optional[ValueWithUnit] = None


class GeometryComputed(BaseModel):
    """Computed geometry values."""
    volume_m3: Optional[float] = None
    surface_m2: Optional[float] = None


class Geometry(BaseModel):
    """Packaging geometry definition."""
    shape: GeometryShape
    dimensions: GeometryDimensions
    computed: Optional[GeometryComputed] = None


# =============================================================================
# Substance Models
# =============================================================================

class SubstanceProperties(BaseModel):
    """Chemical properties of a substance."""
    name: str
    cas: Optional[List[str]] = None
    mw: Optional[float] = Field(None, description="Molecular weight (g/mol)")
    logP: Optional[float] = None
    formula: Optional[str] = None


class RegulatoryStatus(BaseModel):
    """Regulatory status for a region."""
    authorized: Optional[bool] = None
    SML: Optional[float] = Field(None, description="Specific Migration Limit (mg/kg)")
    fcs_listed: Optional[bool] = None
    gb_listed: Optional[bool] = None


class SubstanceRegulatory(BaseModel):
    """Regulatory information across regions."""
    EU: Optional[RegulatoryStatus] = None
    US: Optional[RegulatoryStatus] = None
    CN: Optional[RegulatoryStatus] = None


class Substance(BaseModel):
    """Substance definition."""
    id: str = Field(default_factory=lambda: str(uuid.uuid4())[:8])
    source: SubstanceSource = SubstanceSource.PUBCHEM
    lookup_name: Optional[str] = None  # For PubChem lookup
    cid: Optional[int] = None  # PubChem CID
    properties: Optional[SubstanceProperties] = None
    regulatory: Optional[SubstanceRegulatory] = None
    SML: Optional[float] = Field(None, description="Override SML (mg/kg)")
    color: Optional[str] = Field("#3B82F6", description="Color for plotting")


# =============================================================================
# Layer Models
# =============================================================================

class LayerSubstance(BaseModel):
    """Substance assignment in a layer."""
    substance_id: str
    C0: ValueWithUnit = Field(..., description="Initial concentration")
    D_auto: bool = Field(True, description="Auto-compute diffusivity")
    D_override: Optional[float] = Field(None, description="Manual D override (m2/s)")
    k_auto: bool = Field(True, description="Auto-compute partition coefficient")
    k_override: Optional[float] = Field(None, description="Manual k override")


class Layer(BaseModel):
    """Material layer definition."""
    index: int = Field(..., ge=1, description="Layer index (1 = food contact)")
    polymer: PolymerType
    thickness: ValueWithUnit
    name: Optional[str] = None
    substances: List[LayerSubstance] = Field(default_factory=list)


# =============================================================================
# Food Models
# =============================================================================

class K0Value(BaseModel):
    """Food-side partition coefficient for a substance."""
    value: Optional[float] = None
    auto: bool = True


class Food(BaseModel):
    """Food/contact medium definition."""
    category: FoodCategory = FoodCategory.REALFOOD
    texture: FoodTexture = FoodTexture.LIQUID
    affinity: FoodAffinity = FoodAffinity.FAT
    simulant: Optional[str] = Field(None, description="Simulant type (ethanol, ethanol95, etc.)")
    name: Optional[str] = None
    k0_values: Dict[str, K0Value] = Field(default_factory=dict)


# =============================================================================
# Contact Step Models
# =============================================================================

class ContactStep(BaseModel):
    """Contact/storage step definition."""
    index: int = Field(..., ge=1)
    name: Optional[str] = None
    type: ContactType = ContactType.STORAGE
    temperature: ValueWithUnit
    duration: ValueWithUnit
    with_food_contact: bool = True
    food_override: Optional[Food] = Field(None, description="Override food properties for this step")


# =============================================================================
# Simulation Config Models
# =============================================================================

class OutputUnits(BaseModel):
    """Output units configuration."""
    time: str = "days"
    length: str = "um"
    concentration: str = "mg/kg"


class SimulationConfig(BaseModel):
    """Simulation configuration."""
    solver: str = "senspatankar"
    time_factor: float = Field(2.0, description="Simulate for time_factor * tcontact")
    n_time_points: int = 1000
    n_space_points: int = 200
    output_units: OutputUnits = Field(default_factory=OutputUnits)


# =============================================================================
# Results Models
# =============================================================================

class SubstanceResult(BaseModel):
    """Results for a single substance."""
    substance_id: str
    CF_at_tcontact: float
    CF_equilibrium: float
    SML: Optional[float] = None
    compliant: Optional[bool] = None
    margin_percent: Optional[float] = None


class TimeSeries(BaseModel):
    """Time series data."""
    t_days: List[float]
    CF_by_substance: Dict[str, List[float]]


class ConcentrationProfiles(BaseModel):
    """Concentration profile data."""
    x_um: List[float]
    times_days: List[float]
    Cx_by_substance: Dict[str, List[List[float]]]


class Results(BaseModel):
    """Simulation results."""
    computed_at: Optional[datetime] = None
    elapsed_seconds: Optional[float] = None
    substances: List[SubstanceResult] = Field(default_factory=list)
    time_series: Optional[TimeSeries] = None
    concentration_profiles: Optional[ConcentrationProfiles] = None


# =============================================================================
# Restart State Models
# =============================================================================

class RestartState(BaseModel):
    """State for restarting simulation from previous point."""
    enabled: bool = False
    from_step: int = Field(1, ge=1, description="Contact step to restart from")
    initial_CF: Dict[str, float] = Field(
        default_factory=dict,
        description="Initial CF (mg/kg) per substance"
    )
    initial_Cx: Optional[Dict[str, Dict[str, List[float]]]] = Field(
        None,
        description="Initial Cx profiles per substance per layer"
    )
    source_session: Optional[str] = Field(None, description="Source session file")


# =============================================================================
# Metadata Model
# =============================================================================

class Metadata(BaseModel):
    """Session metadata."""
    name: str
    description: Optional[str] = None
    created: datetime = Field(default_factory=datetime.utcnow)
    modified: datetime = Field(default_factory=datetime.utcnow)
    author: Optional[str] = None
    source_script: Optional[str] = Field(None, description="Original Python script")


# =============================================================================
# Main Session Model
# =============================================================================

class Session(BaseModel):
    """
    Complete SFPPy Studio session.

    This is the root model for .sfppy.json files.
    """
    version: str = "1.0"
    metadata: Metadata
    geometry: Geometry
    substances: List[Substance] = Field(default_factory=list)
    layers: List[Layer] = Field(default_factory=list)
    food: Food = Field(default_factory=Food)
    contact_steps: List[ContactStep] = Field(default_factory=list)
    simulation_config: SimulationConfig = Field(default_factory=SimulationConfig)
    results: Optional[Results] = None
    restart_state: Optional[RestartState] = None

    @field_validator('layers')
    @classmethod
    def validate_layers(cls, v):
        """Ensure at least one layer exists."""
        if not v:
            raise ValueError("At least one layer is required")
        return v

    @field_validator('contact_steps')
    @classmethod
    def validate_contact_steps(cls, v):
        """Ensure at least one contact step exists."""
        if not v:
            raise ValueError("At least one contact step is required")
        return v

    def to_json(self, indent: int = 2) -> str:
        """Serialize to JSON string."""
        return self.model_dump_json(indent=indent, exclude_none=True)

    @classmethod
    def from_json(cls, json_str: str) -> "Session":
        """Deserialize from JSON string."""
        return cls.model_validate_json(json_str)

    def save(self, filepath: str) -> None:
        """Save session to file."""
        from pathlib import Path
        Path(filepath).write_text(self.to_json(), encoding='utf-8')

    @classmethod
    def load(cls, filepath: str) -> "Session":
        """Load session from file."""
        from pathlib import Path
        json_str = Path(filepath).read_text(encoding='utf-8')
        return cls.from_json(json_str)


# =============================================================================
# Batch Processing Model
# =============================================================================

class BatchConfig(BaseModel):
    """Batch processing configuration."""
    name: str
    sessions: List[str] = Field(..., description="List of session file paths")
    parallel: bool = True
    output_folder: str = "results/"


# =============================================================================
# Factory Functions
# =============================================================================

def create_empty_session(name: str = "Untitled Session") -> Session:
    """Create a new empty session with default values."""
    return Session(
        metadata=Metadata(name=name),
        geometry=Geometry(
            shape=GeometryShape.CYLINDER,
            dimensions=GeometryDimensions(
                radius=ValueWithUnit(value=50, unit="mm"),
                height=ValueWithUnit(value=100, unit="mm")
            )
        ),
        layers=[
            Layer(
                index=1,
                polymer=PolymerType.LDPE,
                thickness=ValueWithUnit(value=100, unit="um")
            )
        ],
        contact_steps=[
            ContactStep(
                index=1,
                temperature=ValueWithUnit(value=40, unit="degC"),
                duration=ValueWithUnit(value=10, unit="days")
            )
        ]
    )
