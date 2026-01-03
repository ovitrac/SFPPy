"""
Session Management for SFPPy Studio

Provides in-memory session storage for simulation state including:
- Substances (as migrantToxtree instances)
- Assembly configuration
- Contact conditions
- Computed properties (D, k, C0 matrix)

Author: Olivier Vitrac, PhD, HDR
Organization: INRAE + Generative Simulation
"""

import sys
import uuid
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime, timedelta
from dataclasses import dataclass, field

# Add parent paths
sys.path.insert(0, str(Path(__file__).parent.parent.parent))


# ========== DATA CLASSES ==========

@dataclass
class SubstanceData:
    """Stored substance data from migrantToxtree."""
    id: str
    name: str
    cid: Optional[int]
    cas: Optional[List[str]]  # CAS can be an array
    mw: Optional[float]
    formula: Optional[str]
    smiles: Optional[str]
    logP: Optional[float]

    # Regulatory
    eu_authorized: Optional[bool] = None
    eu_sml: Optional[float] = None
    us_authorized: Optional[bool] = None
    cn_authorized: Optional[bool] = None

    # ToxTree data
    cramer_class: Optional[str] = None  # "I", "II", "III"
    cramer_value: Optional[int] = None  # 0, 1, 2, 3
    ttc: Optional[float] = None  # µg/kg bw/day
    cf_ttc: Optional[float] = None  # mg/kg food intake
    has_alerts: bool = False
    alerts: List[str] = field(default_factory=list)

    # Computed properties (can be updated per temperature)
    D: Optional[float] = None  # Default diffusivity at reference T
    k: Optional[float] = None  # Default partition coefficient

    # Reference to raw migrant object (for property access)
    _migrant_instance: Any = field(default=None, repr=False)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "id": self.id,
            "name": self.name,
            "cid": self.cid,
            "cas": self.cas,
            "mw": self.mw,
            "formula": self.formula,
            "smiles": self.smiles,
            "logP": self.logP,
            "regulatory": {
                "EU": {"authorized": self.eu_authorized, "SML": self.eu_sml},
                "US": {"authorized": self.us_authorized},
                "CN": {"authorized": self.cn_authorized},
            },
            "toxtree": {
                "cramer_class": self.cramer_class,
                "cramer_value": self.cramer_value,
                "ttc": self.ttc,
                "ttc_unit": "µg/kg bw/day",
                "cf_ttc": self.cf_ttc,
                "cf_ttc_unit": "mg/kg food",
                "has_alerts": self.has_alerts,
                "alerts": self.alerts,
            },
            "D": self.D,
            "k": self.k,
        }


@dataclass
class LayerSubstanceAssignment:
    """Assignment of a substance to a layer with concentration."""
    substance_id: str
    layer_index: int
    C0: float  # Initial concentration
    C0_unit: str = "mg/kg"  # mg/kg (ppm) or kg/m³
    D_override: Optional[float] = None
    k_override: Optional[float] = None

    # Computed values (filled when substance properties are computed)
    D_computed: Optional[float] = None
    k_computed: Optional[float] = None
    D_final: Optional[float] = None
    k_final: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "substance_id": self.substance_id,
            "layer_index": self.layer_index,
            "C0": self.C0,
            "C0_unit": self.C0_unit,
            "D_override": self.D_override,
            "k_override": self.k_override,
            "D_computed": self.D_computed,
            "k_computed": self.k_computed,
            "D_final": self.D_final,
            "k_final": self.k_final,
        }


@dataclass
class ContactStep:
    """Contact condition step."""
    index: int
    temperature: float  # Value
    duration: float  # Value
    temperature_unit: str = "degC"  # Unit
    duration_unit: str = "days"  # Unit
    with_food: bool = True  # False for set-off
    simulant: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "index": self.index,
            "temperature": (self.temperature, self.temperature_unit),
            "duration": (self.duration, self.duration_unit),
            "with_food": self.with_food,
            "simulant": self.simulant,
        }


@dataclass
class LayerConfig:
    """Layer configuration."""
    index: int
    polymer: str
    thickness: float
    thickness_unit: str = "um"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "index": self.index,
            "polymer": self.polymer,
            "thickness": (self.thickness, self.thickness_unit),
        }


@dataclass
class SimulationSession:
    """Complete simulation session state."""
    id: str
    created_at: datetime
    updated_at: datetime
    name: str = "Untitled Simulation"

    # Substances in session (keyed by substance_id)
    substances: Dict[str, SubstanceData] = field(default_factory=dict)

    # Layer configuration
    layers: List[LayerConfig] = field(default_factory=list)

    # Substance-to-layer assignments
    assignments: List[LayerSubstanceAssignment] = field(default_factory=list)

    # Contact steps
    contact_steps: List[ContactStep] = field(default_factory=list)

    # Food/contact properties
    food_category: Optional[str] = None
    food_simulant: Optional[str] = None
    food_density: Optional[float] = None  # kg/m³

    # Geometry
    geometry_shape: Optional[str] = None
    geometry_dimensions: Dict[str, float] = field(default_factory=dict)
    volume_m3: Optional[float] = None
    surface_m2: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "name": self.name,
            "created_at": self.created_at.isoformat(),
            "updated_at": self.updated_at.isoformat(),
            "substances": {k: v.to_dict() for k, v in self.substances.items()},
            "substances_count": len(self.substances),
            "layers": [l.to_dict() for l in self.layers],
            "layers_count": len(self.layers),
            "assignments": [a.to_dict() for a in self.assignments],
            "assignments_count": len(self.assignments),
            "contact_steps": [s.to_dict() for s in self.contact_steps],
            "steps_count": len(self.contact_steps),
            "food": {
                "category": self.food_category,
                "simulant": self.food_simulant,
                "density": self.food_density,
            },
            "geometry": {
                "shape": self.geometry_shape,
                "dimensions": self.geometry_dimensions,
                "volume_m3": self.volume_m3,
                "surface_m2": self.surface_m2,
            },
        }

    def is_valid_for_simulation(self) -> tuple:
        """
        Check if session is ready for simulation.

        Returns (is_valid, errors, warnings)
        """
        errors = []
        warnings = []

        # Must have at least 1 layer
        if len(self.layers) < 1:
            errors.append("At least 1 layer is required")

        # Must have at least 1 contact step
        if len(self.contact_steps) < 1:
            errors.append("At least 1 contact step is required")

        # Must have at least 1 substance assigned
        if len(self.assignments) < 1:
            errors.append("At least 1 substance must be assigned to a layer")

        # Check all assigned substances exist
        for assignment in self.assignments:
            if assignment.substance_id not in self.substances:
                errors.append(f"Substance '{assignment.substance_id}' not found in session")

        # Check layer indices in assignments
        layer_indices = {l.index for l in self.layers}
        for assignment in self.assignments:
            if assignment.layer_index not in layer_indices:
                errors.append(f"Layer {assignment.layer_index} not found for substance assignment")

        # Warnings
        if not self.geometry_shape:
            warnings.append("No geometry defined - using default V/S ratio")

        if not self.food_category:
            warnings.append("No food category selected - using default simulant")

        return len(errors) == 0, errors, warnings


# ========== SESSION STORE ==========

class SessionStore:
    """
    In-memory session storage.

    For production, this could be replaced with Redis or database storage.
    Sessions expire after SESSION_TIMEOUT hours of inactivity.
    """

    SESSION_TIMEOUT = timedelta(hours=24)

    def __init__(self):
        self._sessions: Dict[str, SimulationSession] = {}

    def create_session(self, name: str = "Untitled Simulation") -> SimulationSession:
        """Create a new simulation session."""
        session_id = str(uuid.uuid4())
        now = datetime.utcnow()
        session = SimulationSession(
            id=session_id,
            created_at=now,
            updated_at=now,
            name=name,
        )
        self._sessions[session_id] = session
        return session

    def get_session(self, session_id: str) -> Optional[SimulationSession]:
        """Get session by ID, returns None if not found or expired."""
        session = self._sessions.get(session_id)
        if session is None:
            return None

        # Check expiration
        if datetime.utcnow() - session.updated_at > self.SESSION_TIMEOUT:
            del self._sessions[session_id]
            return None

        return session

    def update_session(self, session: SimulationSession) -> SimulationSession:
        """Update session timestamp."""
        session.updated_at = datetime.utcnow()
        self._sessions[session.id] = session
        return session

    def delete_session(self, session_id: str) -> bool:
        """Delete a session."""
        if session_id in self._sessions:
            del self._sessions[session_id]
            return True
        return False

    def list_sessions(self) -> List[Dict[str, Any]]:
        """List all active sessions (summary only)."""
        self._cleanup_expired()
        return [
            {
                "id": s.id,
                "name": s.name,
                "created_at": s.created_at.isoformat(),
                "updated_at": s.updated_at.isoformat(),
                "substances_count": len(s.substances),
                "layers_count": len(s.layers),
            }
            for s in self._sessions.values()
        ]

    def _cleanup_expired(self):
        """Remove expired sessions."""
        now = datetime.utcnow()
        expired = [
            sid for sid, session in self._sessions.items()
            if now - session.updated_at > self.SESSION_TIMEOUT
        ]
        for sid in expired:
            del self._sessions[sid]


# ========== GLOBAL STORE INSTANCE ==========

_session_store = SessionStore()


def get_session_store() -> SessionStore:
    """Get the global session store."""
    return _session_store


# ========== MIGRANT LOADING UTILITIES ==========

def load_migrant_toxtree(query: str) -> Optional[SubstanceData]:
    """
    Load substance data using migrantToxtree.

    Attempts to load via migrant class, then promote to migrantToxtree
    if ToxTree cache is available.

    Args:
        query: Substance name, CAS, or CID

    Returns:
        SubstanceData instance or None if not found
    """
    try:
        from patankar.loadpubchem import migrant, migrantToxtree
    except ImportError:
        return None

    try:
        # Load basic migrant data
        m = migrant(query, verbose=False)
        if m is None or m.cid is None:
            return None

        # Create substance ID
        substance_id = f"cid_{m.cid}" if m.cid else f"name_{query.lower().replace(' ', '_')}"

        # Handle name which can be a list
        name = m.name
        if isinstance(name, list):
            # Prefer the query if it matches
            query_lower = query.lower()
            for n in name:
                if n.lower() == query_lower:
                    name = n
                    break
            else:
                # Use first reasonable name
                name = name[0] if name else query

        # Handle CAS which is a list
        cas_list = None
        if hasattr(m, 'CAS') and m.CAS:
            cas_list = m.CAS if isinstance(m.CAS, list) else [m.CAS]

        # Handle logP which can be array
        logP = None
        if hasattr(m, 'logP') and m.logP is not None:
            logP_val = m.logP
            if hasattr(logP_val, '__iter__') and not isinstance(logP_val, str):
                logP = float(logP_val[0]) if len(logP_val) > 0 else None
            else:
                logP = float(logP_val)

        # Create base substance data
        substance = SubstanceData(
            id=substance_id,
            name=name,
            cid=m.cid,
            cas=cas_list,
            mw=m.M if hasattr(m, 'M') else None,
            formula=m.formula if hasattr(m, 'formula') else None,
            smiles=m.smiles if hasattr(m, 'smiles') else None,
            logP=logP,
            eu_authorized=m.EUauthorized if hasattr(m, 'EUauthorized') else None,
            eu_sml=m.EUSML if hasattr(m, 'EUSML') else None,
            us_authorized=m.USauthorized if hasattr(m, 'USauthorized') else None,
            cn_authorized=m.CNauthorized if hasattr(m, 'CNauthorized') else None,
        )

        # Try to promote to migrantToxtree for ToxTree data
        if hasattr(m, 'ispromovable') and m.ispromovable():
            try:
                mt = m.promote()
                if mt is not None and isinstance(mt, migrantToxtree):
                    # Get ToxTree data
                    if hasattr(mt, 'CramerClass'):
                        substance.cramer_class = mt.CramerClass
                    if hasattr(mt, 'CramerValue'):
                        substance.cramer_value = mt.CramerValue
                    if hasattr(mt, 'TTC') and not isinstance(mt.TTC, list):
                        substance.ttc = mt.TTC
                    if hasattr(mt, 'CFTTC') and not isinstance(mt.CFTTC, list):
                        substance.cf_ttc = mt.CFTTC
                    if hasattr(mt, 'has_alerts'):
                        substance.has_alerts = mt.has_alerts
                    if hasattr(mt, 'showalerts'):
                        alerts = mt.showalerts
                        if isinstance(alerts, dict):
                            substance.alerts = list(alerts.keys())

                    # Store migrant instance for property computation
                    substance._migrant_instance = mt
            except Exception:
                # ToxTree promotion failed, continue with basic data
                substance._migrant_instance = m
        else:
            substance._migrant_instance = m

        return substance

    except Exception:
        return None


def create_unknown_substance(
    name: str = "unknown",
    mw: float = 500.0,
    logP: float = 5.0
) -> SubstanceData:
    """
    Create an 'unknown' substance with default D=1e-12, k=1.

    Used when no specific substance data is available.
    """
    return SubstanceData(
        id=f"unknown_{name.lower().replace(' ', '_')}",
        name=name,
        cid=None,
        cas=None,
        mw=mw,
        formula=None,
        smiles=None,
        logP=logP,
        D=1e-12,  # Default diffusivity m²/s
        k=1.0,    # Default partition coefficient
    )


def compute_substance_properties(
    substance: SubstanceData,
    polymer_code: str,
    temperature_C: float
) -> Dict[str, float]:
    """
    Compute D and k for a substance in a specific polymer at given temperature.

    Uses the SFPPy layer module for proper computation if available.

    Args:
        substance: SubstanceData with migrant instance
        polymer_code: Polymer code (e.g., 'LDPE', 'PP')
        temperature_C: Temperature in Celsius

    Returns:
        Dict with D, k, D0, k0 values
    """
    result = {
        "D": 1e-12,  # Default fallback
        "k": 1.0,
        "D0": None,
        "k0": None,
        "method": "default",
    }

    try:
        from patankar import layer as layer_module

        # Get the polymer class
        polymer_cls = getattr(layer_module, polymer_code, None)
        if polymer_cls is None:
            return result

        # Create layer instance with substance
        migrant_instance = getattr(substance, '_migrant_instance', None)
        if migrant_instance is None:
            return result

        # Create layer with migrant
        try:
            layer_instance = polymer_cls(
                l=(100, "um"),  # Default thickness
                migrant=migrant_instance,
                T=(temperature_C, "degC")
            )

            # Extract computed properties
            # D is the computed diffusivity at temperature T
            if hasattr(layer_instance, 'D') and layer_instance.D is not None:
                D_val = layer_instance.D
                # Handle array values
                if hasattr(D_val, '__iter__') and not isinstance(D_val, str):
                    result["D"] = float(D_val[0]) if len(D_val) > 0 else None
                else:
                    result["D"] = float(D_val)
                result["method"] = "sfppy"

            # D0 is stored in _D (reference/default diffusivity)
            if hasattr(layer_instance, '_D') and layer_instance._D is not None:
                D0_val = layer_instance._D
                if hasattr(D0_val, '__iter__') and not isinstance(D0_val, str):
                    result["D0"] = float(D0_val[0]) if len(D0_val) > 0 else None
                else:
                    result["D0"] = float(D0_val)

            # k is the computed partition coefficient
            if hasattr(layer_instance, 'k') and layer_instance.k is not None:
                k_val = layer_instance.k
                if hasattr(k_val, '__iter__') and not isinstance(k_val, str):
                    result["k"] = float(k_val[0]) if len(k_val) > 0 else None
                else:
                    result["k"] = float(k_val)

            # k0 is stored in _k (reference/default partition coefficient)
            if hasattr(layer_instance, '_k') and layer_instance._k is not None:
                k0_val = layer_instance._k
                if hasattr(k0_val, '__iter__') and not isinstance(k0_val, str):
                    result["k0"] = float(k0_val[0]) if len(k0_val) > 0 else None
                else:
                    result["k0"] = float(k0_val)

        except Exception as e:
            # Fallback to estimation
            result["method"] = "estimation"
            result["error"] = str(e)

    except ImportError:
        result["method"] = "estimation"

    return result


def convert_concentration(
    value: float,
    from_unit: str,
    to_unit: str,
    polymer_density: float = 1000.0  # kg/m³
) -> float:
    """
    Convert concentration between various units.

    Supported units:
    - mg/kg (default, ppm equivalent) - milligrams per kilogram
    - ppm (parts per million = mg/kg)
    - ppb (parts per billion = µg/kg)
    - g/kg (grams per kilogram)
    - kg/kg (mass fraction)
    - ng/kg (nanograms per kilogram)
    - kg/m3 (volumetric, for computation)

    All mass-based units are converted through kg/kg (mass fraction),
    then to kg/m³ using polymer_density for computation.

    Args:
        value: Concentration value
        from_unit: Source unit
        to_unit: Target unit
        polymer_density: Polymer density in kg/m³ (for volumetric conversions)

    Returns:
        Converted value
    """
    if from_unit == to_unit:
        return value

    # Conversion factors TO kg/kg (mass fraction)
    to_kg_kg = {
        "kg/kg": 1.0,
        "g/kg": 1e-3,
        "mg/kg": 1e-6,
        "ppm": 1e-6,      # ppm = mg/kg
        "µg/kg": 1e-9,
        "ug/kg": 1e-9,    # ASCII variant
        "ppb": 1e-9,      # ppb = µg/kg
        "ng/kg": 1e-12,
        "ppt": 1e-12,     # ppt = ng/kg (parts per trillion)
    }

    # Special handling for volumetric units (kg/m³)
    # kg/m³ = (kg/kg) * (kg/m³ density)
    if from_unit == "kg/m3":
        # Convert kg/m³ to kg/kg first
        kg_kg_value = value / polymer_density
    elif from_unit in to_kg_kg:
        kg_kg_value = value * to_kg_kg[from_unit]
    else:
        # Unknown unit, try lowercase
        from_unit_lower = from_unit.lower().replace(" ", "")
        if from_unit_lower in to_kg_kg:
            kg_kg_value = value * to_kg_kg[from_unit_lower]
        else:
            raise ValueError(f"Unknown concentration unit: {from_unit}")

    # Convert from kg/kg to target unit
    if to_unit == "kg/m3":
        # kg/kg * density = kg/m³
        return kg_kg_value * polymer_density
    elif to_unit in to_kg_kg:
        return kg_kg_value / to_kg_kg[to_unit]
    else:
        to_unit_lower = to_unit.lower().replace(" ", "")
        if to_unit_lower in to_kg_kg:
            return kg_kg_value / to_kg_kg[to_unit_lower]
        else:
            raise ValueError(f"Unknown concentration unit: {to_unit}")


# Supported concentration units for UI
CONCENTRATION_UNITS = [
    {"code": "mg/kg", "name": "mg/kg (ppm)", "description": "Milligrams per kilogram", "default": True},
    {"code": "ppm", "name": "ppm", "description": "Parts per million (= mg/kg)"},
    {"code": "g/kg", "name": "g/kg", "description": "Grams per kilogram"},
    {"code": "ppb", "name": "ppb", "description": "Parts per billion (= µg/kg)"},
    {"code": "µg/kg", "name": "µg/kg", "description": "Micrograms per kilogram"},
    {"code": "ng/kg", "name": "ng/kg", "description": "Nanograms per kilogram"},
    {"code": "kg/kg", "name": "kg/kg", "description": "Mass fraction"},
    {"code": "kg/m3", "name": "kg/m³", "description": "Volumetric (requires density)"},
]


def get_concentration_for_computation(
    value: float,
    unit: str,
    polymer_density: float = 1000.0
) -> float:
    """
    Convert concentration to kg/m³ for computation.

    Args:
        value: Concentration value in given unit
        unit: Concentration unit
        polymer_density: Polymer density in kg/m³

    Returns:
        Concentration in kg/m³
    """
    return convert_concentration(value, unit, "kg/m3", polymer_density)


# ========== TTC CLASS VALUES (from migrantToxtree) ==========

TTC_VALUES = [0.0025, 1.5, 9.0, 30]  # µg/kg bw/day for Cramer classes 0, I, II, III
CFTTC_VALUES = [ttc * 60 * 1e-3 for ttc in TTC_VALUES]  # mg/kg food intake
