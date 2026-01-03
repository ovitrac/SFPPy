"""
Simulation Routes - Migration Simulation Engine

API endpoints for running migration simulations with the Patankar solver.
"""

import sys
import uuid
import time
import json
from pathlib import Path
from datetime import datetime
from typing import List, Optional, Dict, Any

from fastapi import APIRouter, HTTPException, BackgroundTasks
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

# Add parent paths
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

router = APIRouter()


# ========== DATA MODELS ==========

class ParameterOverride(BaseModel):
    use_computed: bool = True
    override_value: Optional[float] = None


class LayerConfig(BaseModel):
    index: int = Field(..., ge=1, le=10)
    polymer: str
    thickness: float
    thickness_unit: str = "um"
    D: Optional[float] = Field(None, description="Diffusion coefficient (m²/s)")
    D_override: Optional[ParameterOverride] = None
    k: Optional[float] = Field(None, description="Partition coefficient")
    k_override: Optional[ParameterOverride] = None
    k0: Optional[float] = Field(None, description="Initial partition coefficient")
    k0_override: Optional[ParameterOverride] = None
    C0: float = Field(default=0.0, description="Initial concentration (mg/kg, w:w)")
    # Density for w:w <-> w:v conversion (follows same auto/override pattern as D, k, k0)
    rho: Optional[float] = Field(None, description="Layer density (kg/m³), auto-computed from polymer")
    rho_override: Optional[ParameterOverride] = None


class StepConfig(BaseModel):
    index: int = Field(..., ge=1, le=10)
    temperature_C: float = 25.0
    duration: float
    duration_unit: str = "days"
    with_food: bool = True
    setoff_type: Optional[str] = Field(
        default="stacked",
        description="Set-off configuration when with_food=False: 'stacked' or 'rolled'"
    )


class GeometryConfig(BaseModel):
    shape: str
    dimensions: dict
    volume_m3: Optional[float] = None
    surface_m2: Optional[float] = None


class FoodConfig(BaseModel):
    category: str
    simulant: Optional[str] = None
    CF0: Optional[Dict[str, Any]] = Field(None, description="Initial concentration in food per substance {subId: {value, unit}}")
    # Density for w:w <-> w:v conversion (follows same auto/override pattern as D, k, k0)
    density: Optional[float] = Field(None, description="Food/simulant density (kg/m³), auto-computed from simulant")
    density_override: Optional[ParameterOverride] = None


class SubstanceConfig(BaseModel):
    """Substance configuration for simulation."""
    id: str
    name: str
    cas: Optional[str] = Field(None, description="CAS number for lookup fallback")
    cid: Optional[int] = Field(None, description="PubChem CID for lookup fallback")
    mw: Optional[float] = Field(None, description="Molecular weight (g/mol)")
    logP: Optional[float] = Field(None, description="Log partition coefficient")
    SML: Optional[float] = Field(None, description="Specific Migration Limit (mg/kg food)")
    CF0: Optional[float] = Field(0.0, description="Initial concentration in food (mg/kg)")
    CF0_unit: Optional[str] = Field("mg/kg", description="Unit for CF0")
    layer_assignments: Optional[Dict[int, Dict[str, Any]]] = Field(None, description="Layer-specific C0, D, k values")


class SimulationInput(BaseModel):
    name: str = Field(default="Untitled Simulation")
    layers: List[LayerConfig]
    steps: List[StepConfig]
    geometry: GeometryConfig
    food: FoodConfig
    substances: Optional[List[SubstanceConfig]] = Field(default=None, description="Substances to simulate")
    solver_settings: Optional[Dict[str, Any]] = None


class SimulationResult(BaseModel):
    job_id: str
    status: str
    started_at: str
    completed_at: Optional[str]
    input_config: dict
    results: Optional[dict]
    error: Optional[str]


# ========== SOLVER INTEGRATION ==========

def get_sfppy_modules():
    """
    Get all required SFPPy modules for real simulation.

    Returns dict with:
    - solver: senspatankar function
    - layer: base layer class
    - polymers: dict of polymer classes (LDPE, PP, gPET, etc.)
    - food: food module with all contact classes
    - migrant: migrant class for substance lookup
    - geometry: Packaging3D class
    - setoff_classes: dict mapping setoff type names to classes
    - useroverride: useroverride instance for solver settings
    """
    try:
        from patankar.migration import senspatankar
        from patankar.layer import layer, LDPE, LLDPE, HDPE, PP, gPET, wPET, rPET, PS, PMMA, PA6, PA66, PBT, PEN, HIPS, SBS
        from patankar.loadpubchem import migrant
        from patankar.geometry import Packaging3D
        from patankar.useroverride import useroverride
        import patankar.food as food

        polymers = {
            'LDPE': LDPE, 'LLDPE': LLDPE, 'HDPE': HDPE,
            'PP': PP, 'gPET': gPET, 'wPET': wPET, 'PET': gPET, 'rPET': rPET,
            'PS': PS, 'PMMA': PMMA, 'PA6': PA6, 'PA66': PA66,
            'PBT': PBT, 'PEN': PEN, 'HIPS': HIPS, 'SBS': SBS,
        }

        # Set-off classes for with_food=False steps (stacking/rolling configurations)
        setoff_classes = {
            'stacked': food.stacked,
            'stacking': food.stacked,  # synonym
            'rolled': food.rolled,
            'rolling': food.rolled,    # synonym
            'setoff': food.setoff,     # generic
        }

        return {
            'solver': senspatankar,
            'layer': layer,
            'polymers': polymers,
            'food': food,
            'migrant': migrant,
            'geometry': Packaging3D,
            'setoff_classes': setoff_classes,
            'useroverride': useroverride,
        }
    except ImportError as e:
        print(f"SFPPy import error: {e}")
        return None


def apply_solver_settings(sfppy: dict, solver_settings: Optional[Dict[str, Any]] = None):
    """
    Apply solver settings to useroverride for the Patankar solver.

    Parameters from config/settings that map to useroverride:
    - ntimes: number of stored simulation times (default 1000, max 20000)
    - nmesh: total FV volumes in assembly (default 600)
    - nmeshmin: minimum FV volumes per layer (default 20)
    - timescale: "sqrt" (diffusion-optimal) or "linear"
    - RelTol: relative tolerance (default 1e-6)
    - AbsTol: absolute tolerance (default 1e-6)
    - nmax: number of concentration profiles (default 15)
    """
    useroverride = sfppy.get('useroverride')
    if not useroverride:
        return

    # Load defaults from config if no solver_settings provided
    if solver_settings is None:
        try:
            from app.routes.config import load_config
            config = load_config()
            solver_settings = config.solver.dict()
        except Exception:
            solver_settings = {}

    # Apply settings to useroverride
    if 'ntimes' in solver_settings:
        useroverride.ntimes = min(max(int(solver_settings['ntimes']), 50), 20000)
    if 'nmesh' in solver_settings:
        useroverride.nmesh = min(max(int(solver_settings['nmesh']), 100), 5000)
    if 'nmeshmin' in solver_settings:
        useroverride.nmeshmin = min(max(int(solver_settings['nmeshmin']), 5), 100)
    if 'timescale' in solver_settings:
        useroverride.timescale = solver_settings['timescale'] if solver_settings['timescale'] in ('sqrt', 'linear') else 'sqrt'
    if 'RelTol' in solver_settings:
        useroverride.RelTol = float(solver_settings['RelTol'])
    if 'AbsTol' in solver_settings:
        useroverride.AbsTol = float(solver_settings['AbsTol'])
    if 'nmax' in solver_settings:
        useroverride.nmax = min(max(int(solver_settings['nmax']), 5), 50)


def convert_thickness(value: float, unit: str) -> float:
    """Convert thickness to meters."""
    conversions = {
        "nm": 1e-9, "um": 1e-6, "µm": 1e-6,
        "mm": 1e-3, "cm": 1e-2, "m": 1.0,
    }
    return value * conversions.get(unit.lower(), 1e-6)


def convert_duration(value: float, unit: str) -> float:
    """Convert duration to seconds."""
    conversions = {
        "s": 1, "sec": 1, "seconds": 1,
        "min": 60, "minutes": 60,
        "h": 3600, "hours": 3600,
        "days": 86400, "d": 86400,
        "weeks": 604800, "w": 604800,
        "months": 2592000,
        "years": 31536000,
    }
    return value * conversions.get(unit.lower(), 86400)


# ========== DENSITY CONVERSION HELPERS (w:w <-> w:v) ==========

# Default densities (kg/m³) for common materials when not available from patankar
DEFAULT_DENSITIES = {
    # Polymers
    'LDPE': 920, 'LLDPE': 920, 'HDPE': 960,
    'PP': 910, 'gPET': 1380, 'wPET': 1380, 'PET': 1380, 'rPET': 1380,
    'PS': 1050, 'PMMA': 1190, 'PA6': 1130, 'PA66': 1140,
    'PBT': 1310, 'PEN': 1360, 'HIPS': 1040, 'SBS': 940,
    # Food simulants
    'water': 997, 'ethanol_10': 982, 'ethanol_20': 969, 'ethanol_50': 914,
    'ethanol_95': 789, 'olive_oil': 920, 'iso_octane': 690, 'tenax': 1000,
}


def get_effective_density(
    auto_value: Optional[float],
    override: Optional[ParameterOverride],
    default: float = 1000.0
) -> float:
    """
    Get effective density value using auto/override pattern.

    Priority:
    1. User override if use_computed=False and override_value is set
    2. Auto-computed value from patankar layer/food
    3. Default value (1000 kg/m³ = water density)

    Args:
        auto_value: Density auto-computed from material (kg/m³)
        override: User override settings
        default: Fallback density if nothing else available

    Returns:
        Effective density in kg/m³
    """
    if override and not override.use_computed and override.override_value is not None:
        return override.override_value
    if auto_value is not None:
        return auto_value
    return default


def convert_C0_ww_to_wv(C0_mg_kg: float, rho_kg_m3: float) -> float:
    """
    Convert initial concentration from w:w to w:v.

    w:w (mass/mass): mg migrant per kg polymer matrix
    w:v (mass/volume): mg migrant per L polymer matrix

    Formula: C_wv (mg/L) = C_ww (mg/kg) × ρ (kg/m³) / 1000

    Note: patankar solver expects w:v internally (mg/L or equivalent).
    UI displays w:w (mg/kg) which is more intuitive for users.

    Args:
        C0_mg_kg: Concentration in mg/kg (w:w)
        rho_kg_m3: Material density in kg/m³

    Returns:
        Concentration in mg/L (w:v)

    Example:
        C0 = 1000 mg/kg, ρ = 920 kg/m³ (LDPE)
        C0_wv = 1000 × 920 / 1000 = 920 mg/L
    """
    return C0_mg_kg * rho_kg_m3 / 1000.0


def convert_CF_wv_to_ww(CF_wv: float, rho_kg_m3: float) -> float:
    """
    Convert food concentration from w:v to w:w.

    w:v (mass/volume): mg migrant per L food/simulant
    w:w (mass/mass): mg migrant per kg food/simulant

    Formula: C_ww (mg/kg) = C_wv (mg/L) × 1000 / ρ (kg/m³)

    Note: patankar solver returns CF in w:v (mg/L).
    UI should display CF in w:w (mg/kg) for regulatory compliance.

    Args:
        CF_wv: Concentration in mg/L (w:v)
        rho_kg_m3: Food/simulant density in kg/m³

    Returns:
        Concentration in mg/kg (w:w)

    Example:
        CF = 5 mg/L, ρ = 920 kg/m³ (olive oil)
        CF_ww = 5 × 1000 / 920 = 5.43 mg/kg
    """
    if rho_kg_m3 <= 0:
        rho_kg_m3 = 1000.0  # Fallback to water
    return CF_wv * 1000.0 / rho_kg_m3


def get_layer_density_from_patankar(polymer_name: str, temperature_C: float = 25.0) -> Optional[float]:
    """
    Get layer density from patankar layer instance using density(T) method.

    The patankar layer classes have temperature-dependent density methods, e.g.:
        LDPE.density(T) = 920 * (1 - 3*(T - Td) * 20e-5)

    Args:
        polymer_name: Polymer code (e.g., 'LDPE', 'PP', 'gPET')
        temperature_C: Temperature for T-dependent density (°C)

    Returns:
        Density in kg/m³, or None if not available
    """
    try:
        from patankar.layer import LDPE, LLDPE, HDPE, PP, gPET, wPET, rPET, PS, PMMA, PA6, PA66, PBT, PEN, HIPS, SBS

        polymer_classes = {
            'LDPE': LDPE, 'LLDPE': LLDPE, 'HDPE': HDPE,
            'PP': PP, 'gPET': gPET, 'wPET': wPET, 'PET': gPET, 'rPET': rPET,
            'PS': PS, 'PMMA': PMMA, 'PA6': PA6, 'PA66': PA66,
            'PBT': PBT, 'PEN': PEN, 'HIPS': HIPS, 'SBS': SBS,
        }

        PolymerClass = polymer_classes.get(polymer_name)
        if PolymerClass:
            # Create instance with temperature to access density(T) method
            layer_instance = PolymerClass(l=100e-6, T=(temperature_C, 'degC'))
            # Use density(T) method which incorporates temperature dependence
            if hasattr(layer_instance, 'density') and callable(layer_instance.density):
                rho_result = layer_instance.density(temperature_C)
                # density() returns tuple (array([920.]), "kg/m**3")
                if isinstance(rho_result, tuple):
                    rho_val = rho_result[0]
                    # Handle numpy array: array([920.]) -> 920.0
                    if hasattr(rho_val, '__iter__') and not isinstance(rho_val, str):
                        return float(rho_val[0])
                    return float(rho_val)
                return float(rho_result)
            # Fallback to rho attribute if density() method not available
            elif hasattr(layer_instance, 'rho'):
                rho = layer_instance.rho
                if hasattr(rho, '__iter__') and not isinstance(rho, str):
                    return float(rho[0])
                return float(rho)
    except Exception as e:
        print(f"Warning: Could not get density for {polymer_name}: {e}")

    # Fallback to default densities
    return DEFAULT_DENSITIES.get(polymer_name)


def get_food_density_from_simulant(simulant: str, temperature_C: float = 25.0) -> float:
    """
    Get food/simulant density.

    Args:
        simulant: Simulant code (e.g., 'ethanol_50', 'olive_oil', 'water')
        temperature_C: Temperature for T-dependent density

    Returns:
        Density in kg/m³ (defaults to 1000 if not found)
    """
    # Map common simulant names to density
    simulant_densities = {
        'water': 997.0,
        'ethanol': 789.0,  # Pure ethanol
        'ethanol_10': 982.0,
        'ethanol_20': 969.0,
        'ethanol_50': 914.0,
        'ethanol_95': 789.0,
        'olive_oil': 920.0,
        'vegetable_oil': 920.0,
        'sunflower_oil': 920.0,
        'iso_octane': 690.0,
        'isooctane': 690.0,
        'tenax': 1000.0,  # Solid adsorbent, treated as unit density
        'mppo': 1060.0,  # Modified PPO
    }

    # Normalize simulant name
    simulant_lower = simulant.lower().replace('-', '_').replace(' ', '_')

    # Direct match
    if simulant_lower in simulant_densities:
        rho_25 = simulant_densities[simulant_lower]
    else:
        # Try partial match
        for key, rho in simulant_densities.items():
            if key in simulant_lower or simulant_lower in key:
                rho_25 = rho
                break
        else:
            rho_25 = 1000.0  # Default to water

    # Simple T-dependence: ρ(T) ≈ ρ_25 × (1 - α × (T - 25))
    # α ≈ 0.001 /°C for most liquids
    alpha = 0.001
    rho_T = rho_25 * (1 - alpha * (temperature_C - 25))

    return max(rho_T, 500.0)  # Sanity bound


# ========== SIMULATION STORAGE ==========

JOBS_DIR = Path(__file__).parent.parent.parent / "jobs"
JOBS_DIR.mkdir(exist_ok=True)

# In-memory cache for active simulations
ACTIVE_SIMULATIONS: Dict[str, Dict[str, Any]] = {}


def save_job(job_id: str, data: dict):
    """Save job data to disk."""
    job_file = JOBS_DIR / f"{job_id}.json"
    with open(job_file, 'w') as f:
        json.dump(data, f, indent=2, default=str)


def load_job(job_id: str) -> Optional[dict]:
    """Load job data from disk."""
    job_file = JOBS_DIR / f"{job_id}.json"
    if job_file.exists():
        with open(job_file, 'r') as f:
            return json.load(f)
    return None


# ========== SIMULATION ENGINE ==========

def run_simulation_task(job_id: str, config: SimulationInput):
    """
    Run simulation as a background task.

    This function interfaces with the SFPPy Patankar solver.
    """
    start_time = time.time()

    try:
        # Update status
        ACTIVE_SIMULATIONS[job_id]["status"] = "running"

        sfppy = get_sfppy_modules()

        if sfppy is None:
            # Fallback to mock results for testing
            results = generate_mock_results(config)
        else:
            # Real simulation using Patankar solver
            results = execute_real_simulation(config, sfppy)

        # Store results
        elapsed = time.time() - start_time
        ACTIVE_SIMULATIONS[job_id].update({
            "status": "completed",
            "completed_at": datetime.utcnow().isoformat() + "Z",
            "results": results,
            "elapsed_seconds": elapsed,
        })

        # Save to disk
        save_job(job_id, ACTIVE_SIMULATIONS[job_id])

    except Exception as e:
        import traceback
        ACTIVE_SIMULATIONS[job_id].update({
            "status": "failed",
            "completed_at": datetime.utcnow().isoformat() + "Z",
            "error": str(e),
            "traceback": traceback.format_exc(),
        })
        save_job(job_id, ACTIVE_SIMULATIONS[job_id])


def execute_real_simulation(config: SimulationInput, sfppy: dict) -> dict:
    """
    Execute real simulation using SFPPy Patankar solver with proper multi-step handling.

    This function implements chained simulations following Example 3 pattern:
    - Each step has its own temperature (D(T) is recalculated via Piringer/Welle models)
    - Steps with with_food=False use set-off boundary conditions (stacked/rolled)
    - C(x) from step n becomes C0 for step n+1
    - Kinetics are combined across all steps

    Solver parameters are applied via useroverride:
    - ntimes, nmesh, nmeshmin: resolution control
    - timescale: "sqrt" for diffusion-optimal time sampling
    - RelTol, AbsTol: ODE integration tolerances

    Args:
        config: SimulationInput with layers, steps, geometry, food, substances
        sfppy: dict from get_sfppy_modules() with solver, polymers, food, etc.

    Returns:
        dict with combined kinetics, per-step profiles, and compliance assessment
    """
    import numpy as np

    # Apply solver settings from config or input
    apply_solver_settings(sfppy, config.solver_settings)

    solver = sfppy['solver']
    polymers = sfppy['polymers']
    food_module = sfppy['food']
    migrant_cls = sfppy['migrant']
    Packaging3D = sfppy['geometry']
    setoff_classes = sfppy['setoff_classes']

    # Build geometry
    geom = config.geometry
    dims = geom.dimensions
    shape = geom.shape.lower()

    if shape == 'cylinder':
        packaging = Packaging3D('Cylinder',
            radius=(dims.get('radius', 50), 'mm'),
            length=(dims.get('height', 100), 'mm'))
    elif shape == 'bottle':
        packaging = Packaging3D('bottle',
            body_radius=(dims.get('body_radius', 40), 'mm'),
            body_height=(dims.get('body_height', 200), 'mm'),
            neck_radius=(dims.get('neck_radius', 15), 'mm'),
            neck_height=(dims.get('neck_height', 50), 'mm'))
    elif shape in ('box', 'rectangle', 'box_container'):
        packaging = Packaging3D('box_container',
            length=(dims.get('length', 200), 'mm'),
            width=(dims.get('width', 150), 'mm'),
            height=(dims.get('height', 50), 'mm'))
    else:
        packaging = Packaging3D('Cylinder', radius=(50, 'mm'), length=(100, 'mm'))

    volume, surface = packaging.get_volume_and_area()

    # Build food class based on category (used for with_food=True steps)
    food_category = config.food.category.lower()
    simulant = config.food.simulant or 'ethanol'

    food_bases = [food_module.realfood]
    if 'liquid' in food_category:
        food_bases.append(food_module.liquid)
    else:
        food_bases.append(food_module.semisolid)

    if 'fat' in food_category or food_category == 'fatty':
        food_bases.append(food_module.fat)
    elif 'aqueous' in food_category:
        food_bases.append(food_module.aqueous)
    else:
        food_bases.append(food_module.fat)

    FoodClass = type('DynamicFood', tuple(food_bases), {'name': 'food'})

    # Get food/simulant density for output conversion (w:v → w:w)
    # Use auto/override pattern
    food_density_auto = get_food_density_from_simulant(simulant, temperature_C=25.0)
    food_density = get_effective_density(
        config.food.density,
        config.food.density_override,
        default=food_density_auto
    )

    # Sort steps by index
    sorted_steps = sorted(config.steps, key=lambda s: s.index)

    # Calculate total contact time (all steps) for result structure
    total_contact_time_s = sum(
        convert_duration(s.duration, s.duration_unit)
        for s in sorted_steps
    )

    # Colors for multi-substance plotting
    substance_colors = [
        '#3B82F6', '#10B981', '#F59E0B', '#EF4444',
        '#8B5CF6', '#EC4899', '#06B6D4', '#84CC16',
    ]

    # Get substances to simulate
    substances_config = config.substances or []

    # If no substances defined, create a default one
    if not substances_config:
        # Get C0 from first layer for default substance
        default_C0 = config.layers[0].C0 if config.layers and config.layers[0].C0 else 1000
        substances_config = [SubstanceConfig(
            id="default",
            name="Substance",
            SML=6.0,
            CF0=0,
            layer_assignments={config.layers[0].index: {'C0': default_C0}} if config.layers else None
        )]

    results_by_substance = []
    overall_compliant = True
    all_x = None
    warnings = []

    # CAS number mapping for common substances (fallback for name lookup)
    CAS_ALIASES = {
        'irganox 1076': '2082-79-3',
        'irganox 1010': '6683-19-8',
        'bht': '128-37-0',
        'butylated hydroxytoluene': '128-37-0',
        'deha': '103-23-1',
        'dehp': '117-81-7',
        'uvitex ob': '7128-64-5',
        'tinuvin 326': '3896-11-5',
        'erucamide': '112-84-5',
    }

    # Process each substance
    for sub_idx, sub_config in enumerate(substances_config):
        # Get migrant - CAS-RN is the primary lookup key (most reliable)
        # Lookup order: CAS → name → CID → known aliases
        m = None
        lookup_attempts = []

        # Try 1: CAS number lookup (PRIMARY KEY - most reliable)
        if sub_config.cas:
            try:
                m = migrant_cls(sub_config.cas)
                # migrant() may return None if not found, or raise ValueError
            except Exception:
                m = None
            lookup_attempts.append(f"CAS='{sub_config.cas}'")

        # Try 2: Direct name lookup (brand names, acronyms, IUPAC names)
        if m is None:
            try:
                m = migrant_cls(sub_config.name)
            except Exception:
                m = None
            lookup_attempts.append(f"name='{sub_config.name}'")

        # Try 3: CID lookup if name failed
        if m is None and sub_config.cid:
            try:
                m = migrant_cls(str(sub_config.cid))
            except Exception:
                m = None
            lookup_attempts.append(f"CID={sub_config.cid}")

        # Try 4: Use CAS alias for common substances (fallback for known additives)
        if m is None:
            name_lower = sub_config.name.lower().strip()
            cas_alias = CAS_ALIASES.get(name_lower)
            if cas_alias:
                try:
                    m = migrant_cls(cas_alias)
                except Exception:
                    m = None
                lookup_attempts.append(f"CAS_alias='{cas_alias}'")

        # If still no migrant, provide clear error message
        if m is None:
            error_msg = (
                f"Substance '{sub_config.name}' not found in PubChem. "
                f"Tried: {', '.join(lookup_attempts)}. "
                f"Please verify the CAS number, name (IUPAC, brand, or acronym), or CID is correct."
            )
            warnings.append(error_msg)
            continue  # Skip to next substance

        # Double-check migrant is valid before proceeding
        if not hasattr(m, 'cid'):
            warnings.append(f"Migrant for '{sub_config.name}' is invalid (no CID attribute) - skipping")
            continue

        # Build initial multilayer structure ONCE
        # The >> operator will modify C0 internally between steps
        layers_list = []
        layer_densities = {}  # Store layer densities for reference
        for lc in sorted(config.layers, key=lambda x: x.index):
            polymer_name = lc.polymer if lc.polymer else 'LDPE'
            PolymerClass = polymers.get(polymer_name, polymers['LDPE'])
            thickness = convert_thickness(lc.thickness, lc.thickness_unit)

            # Get layer density using auto/override pattern (stored for traceability)
            layer_density_auto = get_layer_density_from_patankar(polymer_name)
            layer_density = get_effective_density(
                lc.rho if lc.rho is not None else layer_density_auto,
                lc.rho_override,
                default=DEFAULT_DENSITIES.get(polymer_name, 920.0)
            )
            layer_densities[lc.index] = layer_density

            # Get initial C0 for this substance in this layer (in mg/kg, w:w)
            # Note: patankar solver treats C0 as "arbitrary units" (see layer.py)
            # The convention in SFPPy examples is to use mg/kg (w:w) consistently
            # No conversion is needed - solver is unit-agnostic
            C0 = 0
            if sub_config.layer_assignments and lc.index in sub_config.layer_assignments:
                C0 = sub_config.layer_assignments[lc.index].get('C0', 0)

            layer_obj = PolymerClass(
                l=thickness,
                C0=C0,  # mg/kg (w:w) - solver treats as arbitrary units
                substance=m  # Guaranteed non-None since we continue above if m is None
            )
            layers_list.append(layer_obj)

        # Combine layers (first layer is contact layer)
        multilayer = layers_list[0]
        for l in layers_list[1:]:
            multilayer = multilayer + l

        # Get initial CF0 for this substance
        cf0_value = sub_config.CF0 or 0
        if cf0_value == 0 and config.food.CF0:
            cf0_data = config.food.CF0.get(sub_config.id, {})
            cf0_value = cf0_data.get('value', 0) if isinstance(cf0_data, dict) else (cf0_data or 0)

        # Chain simulations through all steps using >> operator pattern
        # This mirrors the SFPPy core: kin = F.migration(layer), then kin2 = kin >> F2
        step_results = []
        step_profiles = []
        step_simulations = []  # Store simulation objects for combining
        cumulative_time_s = 0
        previous_result = None  # Track previous simulation result for chaining

        for step_idx, step in enumerate(sorted_steps):
            step_temp = (step.temperature_C, 'degC')
            step_duration_s = convert_duration(step.duration, step.duration_unit)

            if step.with_food:
                # Create food contact condition
                medium = FoodClass(
                    volume=volume,
                    surfacearea=surface,
                    contacttime=(step_duration_s, 's'),
                    contacttemperature=step_temp,
                    simulant=simulant,
                    substance=m,  # Guaranteed non-None since we continue above if m is None
                    CF0=cf0_value if step_idx == 0 else 0  # CF0 only for first food contact
                )
                step_type = "food_contact"
            else:
                # Create set-off condition (stacked/rolled)
                setoff_type = getattr(step, 'setoff_type', 'stacked') or 'stacked'
                SetoffClass = setoff_classes.get(setoff_type.lower(), setoff_classes['stacked'])
                medium = SetoffClass(
                    contacttime=(step_duration_s, 's'),
                    contacttemperature=step_temp
                )
                step_type = f"setoff_{setoff_type}"

            # Propagate geometry to medium (packaging >> medium)
            packaging >> medium

            # Chain simulations properly:
            # - Step 1: medium >> multilayer >> medium (runs initial simulation)
            # - Step 2+: previous_result >> medium (chains from previous state with CF0 propagation)
            try:
                if previous_result is None:
                    # First step: run initial simulation
                    medium >> multilayer  # Temperature transfer
                    multilayer >> medium  # Run simulation via >> operator
                    result = medium.lastsimulation  # Get result from lastsimulation
                else:
                    # Subsequent steps: chain from previous result
                    # This calls previous_result.__rshift__(medium) which:
                    # 1. Transfers temperature via medium >> self._lastmultilayer
                    # 2. Calls chaining() with CF0=self.restart.CF to continue simulation
                    result = previous_result >> medium

                # Update previous_result for next iteration
                previous_result = result

                # Store simulation object for combining via + operator
                step_simulations.append(result)

                # Get CF value (0 for setoff steps)
                # Use CFtarget (concentration at contact time) rather than CF[-1] (last array element)
                # Note: patankar solver treats concentrations as "arbitrary units" (see layer.py)
                # SFPPy convention uses mg/kg (w:w) throughout - no conversion needed
                step_CF = result.CF.tolist() if hasattr(result, 'CF') and result.CF is not None else []

                if hasattr(result, 'CFtarget') and result.CFtarget is not None:
                    step_CF_final = float(result.CFtarget)
                else:
                    step_CF_final = float(result.CF[-1]) if step_CF else 0

                # Store step result
                step_results.append({
                    "step_index": step.index,
                    "step_type": step_type,
                    "temperature_C": step.temperature_C,
                    "duration_s": step_duration_s,
                    "duration_days": step_duration_s / 86400,
                    "with_food": step.with_food,
                    "t_s": (result.t + cumulative_time_s).tolist(),
                    "t_days": ((result.t + cumulative_time_s) / 86400).tolist(),
                    "CF_mg_kg": step_CF,  # Now in w:w (mg/kg)
                    "CF_final": step_CF_final,  # Now in w:w (mg/kg)
                })

                # Extract concentration profiles for this step
                if hasattr(result, 'Cx') and result.Cx is not None:
                    n_profiles = 5  # Fewer profiles per step for clarity
                    t_max = result.t[-1]
                    sqrt_t_targets = np.linspace(0, np.sqrt(t_max), n_profiles)
                    t_targets = sqrt_t_targets ** 2

                    profile_indices = [min(np.searchsorted(result.t, t_val), len(result.t) - 1) for t_val in t_targets]
                    step_profile_times = [(result.t[idx] + cumulative_time_s) / 86400 for idx in profile_indices]
                    step_profile_Cx = [result.Cx[idx, :].tolist() for idx in profile_indices]

                    step_profiles.append({
                        "step_index": step.index,
                        "step_type": step_type,
                        "temperature_C": step.temperature_C,
                        "times_days": step_profile_times,
                        "Cx_mg_kg": step_profile_Cx,
                    })

                # Store x coordinates
                if all_x is None and hasattr(result, 'x'):
                    all_x = result.x

                cumulative_time_s += step_duration_s

            except Exception as e:
                import traceback
                warnings.append(f"Step {step.index} simulation failed: {str(e)}")
                traceback.print_exc()
                cumulative_time_s += step_duration_s
                continue

        # Combine step simulations using + operator (like: sol123 = medium1.lastsimulation + medium2.lastsimulation + medium3.lastsimulation)
        combined_sim = None
        if step_simulations:
            combined_sim = step_simulations[0]
            for sim in step_simulations[1:]:
                combined_sim = combined_sim + sim

        # Use combined simulation for unified kinetics (proper time alignment and CF tracking)
        # Note: patankar solver treats concentrations as "arbitrary units" (see layer.py)
        # SFPPy convention uses mg/kg (w:w) throughout - no conversion needed
        if combined_sim is not None:
            combined_t_s = combined_sim.t.tolist()
            combined_t_days = (combined_sim.t / 86400).tolist()
            combined_CF = combined_sim.CF.tolist() if hasattr(combined_sim, 'CF') and combined_sim.CF is not None else []
            # Use CFtarget for correct concentration at contact time
            if hasattr(combined_sim, 'CFtarget') and combined_sim.CFtarget is not None:
                final_CF = float(combined_sim.CFtarget)
            else:
                final_CF = float(combined_sim.CF[-1]) if combined_CF else 0
        else:
            # Fallback to manual combining if + operator fails
            combined_t_s = []
            combined_t_days = []
            combined_CF = []

            for sr in step_results:
                combined_t_s.extend(sr["t_s"])
                combined_t_days.extend(sr["t_days"])
                if sr["CF_mg_kg"]:
                    combined_CF.extend(sr["CF_mg_kg"])

            final_CF = combined_CF[-1] if combined_CF else 0
        SML = sub_config.SML or 6.0
        compliant = final_CF < SML

        # Combine concentration profiles from all steps
        combined_profile_times = []
        combined_profile_Cx = []
        for sp in step_profiles:
            combined_profile_times.extend(sp["times_days"])
            combined_profile_Cx.extend(sp["Cx_mg_kg"])

        results_by_substance.append({
            "id": sub_config.id,
            "name": sub_config.name,
            "color": substance_colors[sub_idx % len(substance_colors)],
            "CF_mg_kg": combined_CF,  # w:w (mg/kg)
            "CF_at_tcontact": final_CF,  # w:w (mg/kg)
            "equilibrium_CF_mg_kg": final_CF,  # w:w (mg/kg)
            "SML_mg_kg": SML,
            "compliant": compliant,
            "margin_percent": ((SML - final_CF) / SML * 100) if compliant else 0,
            "steps": step_results,
            "concentration_profile": {
                "times_days": combined_profile_times,
                "Cx_mg_kg": combined_profile_Cx,
            } if combined_profile_Cx else None,
            "step_profiles": step_profiles,
            # Density information for traceability
            "layer_densities_kg_m3": layer_densities,
            "food_density_kg_m3": food_density,
        })

        if not compliant:
            overall_compliant = False

    # Build response
    x_um = (all_x * 1e6).tolist() if all_x is not None else []
    first_sub = results_by_substance[0] if results_by_substance else {}

    return {
        "time_s": first_sub.get("CF_mg_kg", []),  # Placeholder - actual times
        "time_days": [t / 86400 for t in range(0, int(total_contact_time_s), max(1, int(total_contact_time_s / 200)))],
        "tcontact_s": total_contact_time_s,
        "tcontact_days": total_contact_time_s / 86400,
        "substances": results_by_substance,
        "CF_mg_kg": first_sub.get("CF_mg_kg", []),
        "final_CF_mg_kg": first_sub.get("CF_at_tcontact", 0),
        "equilibrium_CF_mg_kg": first_sub.get("equilibrium_CF_mg_kg", 0),
        "SML_mg_kg": first_sub.get("SML_mg_kg"),
        "compliant": overall_compliant,
        "margin_percent": first_sub.get("margin_percent"),
        "concentration_profiles": {
            "x_um": x_um,
            "times_days": first_sub.get("concentration_profile", {}).get("times_days", []) if first_sub.get("concentration_profile") else [],
            "Cx_mg_kg": first_sub.get("concentration_profile", {}).get("Cx_mg_kg", []) if first_sub.get("concentration_profile") else [],
        },
        "steps_summary": [
            {
                "index": s.index,
                "temperature_C": s.temperature_C,
                "duration": s.duration,
                "duration_unit": s.duration_unit,
                "with_food": s.with_food,
                "setoff_type": getattr(s, 'setoff_type', 'stacked') if not s.with_food else None,
            }
            for s in sorted_steps
        ],
        # Density information (conversion NOT applied - solver uses arbitrary units)
        # SFPPy convention uses mg/kg (w:w) consistently for C0 and CF
        "density_info": {
            "conversion_applied": False,
            "unit": "mg/kg (w:w)",
            "note": "Patankar solver treats concentrations as arbitrary units; SFPPy convention is mg/kg",
            "food_density_kg_m3": food_density,
            "simulant": simulant,
        },
        "warnings": warnings,
        "is_mock": False,
        "multi_step": len(sorted_steps) > 1,
    }


def generate_mock_results(config: SimulationInput) -> dict:
    """
    Generate mock results for testing when SFPPy is not available.

    Following migration.py conventions:
    - Simulate for 2*tcontact to ensure data beyond target
    - Report values at tcontact
    - Use sqrt time scale for dense early-time sampling
    - Support multiple substances with individual SML values
    """
    import math
    import numpy as np

    # Calculate contact time (sum of all steps with food contact)
    tcontact_s = sum(
        convert_duration(s.duration, s.duration_unit)
        for s in config.steps if s.with_food
    )

    # Simulate for 2*tcontact (standard in migration.py)
    tmax_s = 2 * tcontact_s

    # Generate time points using sqrt spacing (denser at early times)
    n_points = 200
    sqrt_t = np.linspace(0, np.sqrt(tmax_s), n_points)
    times = sqrt_t ** 2
    times = times.tolist()

    # Get substances from config or use default
    substances = config.substances or []

    # Colors for multi-substance plotting (consistent with migration.py colormap)
    substance_colors = [
        '#3B82F6',  # blue
        '#10B981',  # emerald
        '#F59E0B',  # amber
        '#EF4444',  # red
        '#8B5CF6',  # violet
        '#EC4899',  # pink
        '#06B6D4',  # cyan
        '#84CC16',  # lime
    ]

    results_by_substance = []
    overall_compliant = True

    if not substances:
        # Single default substance (backward compatibility)
        # Get C0 from first layer
        C0 = config.layers[0].C0 if config.layers and config.layers[0].C0 else 1000

        # Simple diffusion model: CF(t) = CF_eq * (1 - exp(-t/tau))
        # CF_eq depends on C0, layer thickness, and V/S ratio
        CF_eq = C0 * 0.5  # mg/kg in food (simplified)
        tau = tcontact_s / 3  # time constant

        CF_values = [CF_eq * (1 - math.exp(-t / tau)) if tau > 0 else 0 for t in times]

        # Find CF at tcontact (not at 2*tcontact)
        tcontact_idx = np.searchsorted(times, tcontact_s)
        if tcontact_idx >= len(CF_values):
            tcontact_idx = len(CF_values) - 1
        CF_at_tcontact = CF_values[tcontact_idx]

        SML = 6.0  # mg/kg food (default)
        compliant = CF_at_tcontact < SML if SML else True

        results_by_substance.append({
            "id": "default",
            "name": "Substance",
            "color": substance_colors[0],
            "CF_mg_kg": CF_values,
            "CF_at_tcontact": CF_at_tcontact,
            "equilibrium_CF_mg_kg": CF_eq,
            "SML_mg_kg": SML,
            "compliant": compliant,
            "margin_percent": ((SML - CF_at_tcontact) / SML * 100) if SML and SML > 0 else None,
        })
        overall_compliant = compliant
    else:
        # Multiple substances
        for i, sub in enumerate(substances):
            # Get C0 for this substance (from layer_assignments or default)
            C0 = 1000  # default
            if sub.layer_assignments:
                # Sum C0 from all layers
                C0 = sum(
                    la.get('C0', 0)
                    for la in sub.layer_assignments.values()
                )

            # Vary tau slightly per substance (based on molecular weight if available)
            tau_factor = 1 + (i * 0.1)  # slight variation for visual distinction
            tau = tcontact_s / (3 * tau_factor)

            CF_eq = C0 * (0.4 + i * 0.05)  # varied equilibrium
            CF_values = [CF_eq * (1 - math.exp(-t / tau)) if tau > 0 else 0 for t in times]

            # Find CF at tcontact
            tcontact_idx = np.searchsorted(times, tcontact_s)
            if tcontact_idx >= len(CF_values):
                tcontact_idx = len(CF_values) - 1
            CF_at_tcontact = CF_values[tcontact_idx]

            SML = sub.SML or 6.0
            compliant = CF_at_tcontact < SML if SML else True

            color = substance_colors[i % len(substance_colors)]

            results_by_substance.append({
                "id": sub.id,
                "name": sub.name,
                "color": color,
                "CF_mg_kg": CF_values,
                "CF_at_tcontact": CF_at_tcontact,
                "equilibrium_CF_mg_kg": CF_eq,
                "SML_mg_kg": SML,
                "compliant": compliant,
                "margin_percent": ((SML - CF_at_tcontact) / SML * 100) if SML and SML > 0 else None,
            })

            if not compliant:
                overall_compliant = False

    # Generate mock concentration profiles (Cx vs x at selected times)
    n_x = 50
    x_positions = np.linspace(0, 1, n_x).tolist()  # normalized 0-1

    # Select sqrt-spaced time indices for profiles
    n_profiles = 10
    profile_times = []
    profile_Cx = []

    sqrt_profile_t = np.linspace(0, np.sqrt(tcontact_s), n_profiles)
    profile_t_values = (sqrt_profile_t ** 2).tolist()

    for t_val in profile_t_values:
        profile_times.append(t_val / 86400)  # convert to days
        # Mock concentration profile: exponential decay from layer boundary
        decay_factor = 1 - math.exp(-t_val / (tcontact_s / 2)) if tcontact_s > 0 else 0
        Cx = [1000 * (1 - decay_factor * (1 - math.exp(-5 * x))) for x in x_positions]
        profile_Cx.append(Cx)

    return {
        "time_s": times,
        "time_days": [t / 86400 for t in times],
        "tcontact_s": tcontact_s,
        "tcontact_days": tcontact_s / 86400,
        "substances": results_by_substance,
        # Backward compatibility: first substance values
        "CF_mg_kg": results_by_substance[0]["CF_mg_kg"] if results_by_substance else [],
        "final_CF_mg_kg": results_by_substance[0]["CF_at_tcontact"] if results_by_substance else 0,
        "equilibrium_CF_mg_kg": results_by_substance[0]["equilibrium_CF_mg_kg"] if results_by_substance else 0,
        "SML_mg_kg": results_by_substance[0]["SML_mg_kg"] if results_by_substance else None,
        "compliant": overall_compliant,
        "margin_percent": results_by_substance[0]["margin_percent"] if results_by_substance else None,
        # Concentration profiles
        "concentration_profiles": {
            "x_normalized": x_positions,
            "times_days": profile_times,
            "Cx_mg_kg": profile_Cx,
        },
        "warnings": [],
        "is_mock": True,
    }


# ========== API ENDPOINTS ==========

@router.post("/run")
async def run_simulation(
    config: SimulationInput,
    background_tasks: BackgroundTasks,
    async_mode: bool = True,
):
    """
    Start a migration simulation.

    Args:
        config: Complete simulation configuration
        async_mode: If True, run in background and return job_id immediately
    """
    job_id = str(uuid.uuid4())[:8]

    # Initialize job
    job_data = {
        "job_id": job_id,
        "name": config.name,
        "status": "pending",
        "created_at": datetime.utcnow().isoformat() + "Z",
        "started_at": None,
        "completed_at": None,
        "input_config": config.dict(),
        "results": None,
        "error": None,
    }

    ACTIVE_SIMULATIONS[job_id] = job_data
    save_job(job_id, job_data)

    if async_mode:
        # Run in background
        job_data["status"] = "queued"
        job_data["started_at"] = datetime.utcnow().isoformat() + "Z"
        background_tasks.add_task(run_simulation_task, job_id, config)

        return JSONResponse({
            "success": True,
            "job_id": job_id,
            "status": "queued",
            "message": "Simulation queued for execution",
        })
    else:
        # Run synchronously
        job_data["started_at"] = datetime.utcnow().isoformat() + "Z"
        run_simulation_task(job_id, config)

        return JSONResponse({
            "success": ACTIVE_SIMULATIONS[job_id]["status"] == "completed",
            "job_id": job_id,
            "status": ACTIVE_SIMULATIONS[job_id]["status"],
            "results": ACTIVE_SIMULATIONS[job_id].get("results"),
            "error": ACTIVE_SIMULATIONS[job_id].get("error"),
        })


@router.get("/status/{job_id}")
async def get_simulation_status(job_id: str):
    """Get status of a simulation job."""
    # Check in-memory first
    if job_id in ACTIVE_SIMULATIONS:
        job = ACTIVE_SIMULATIONS[job_id]
        return JSONResponse({
            "success": True,
            "job_id": job_id,
            "status": job["status"],
            "created_at": job.get("created_at"),
            "started_at": job.get("started_at"),
            "completed_at": job.get("completed_at"),
            "has_results": job.get("results") is not None,
        })

    # Check disk
    job = load_job(job_id)
    if job:
        ACTIVE_SIMULATIONS[job_id] = job  # Cache it
        return JSONResponse({
            "success": True,
            "job_id": job_id,
            "status": job["status"],
            "created_at": job.get("created_at"),
            "started_at": job.get("started_at"),
            "completed_at": job.get("completed_at"),
            "has_results": job.get("results") is not None,
        })

    raise HTTPException(status_code=404, detail=f"Job {job_id} not found")


@router.get("/results/{job_id}")
async def get_simulation_results(job_id: str):
    """Get full results of a completed simulation."""
    # Check in-memory first
    if job_id in ACTIVE_SIMULATIONS:
        job = ACTIVE_SIMULATIONS[job_id]
    else:
        job = load_job(job_id)
        if job:
            ACTIVE_SIMULATIONS[job_id] = job

    if not job:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    if job["status"] not in ["completed", "failed"]:
        return JSONResponse({
            "success": False,
            "job_id": job_id,
            "status": job["status"],
            "message": "Simulation not yet completed",
        })

    return JSONResponse({
        "success": job["status"] == "completed",
        "job_id": job_id,
        "status": job["status"],
        "name": job.get("name"),
        "created_at": job.get("created_at"),
        "completed_at": job.get("completed_at"),
        "elapsed_seconds": job.get("elapsed_seconds"),
        "results": job.get("results"),
        "error": job.get("error"),
    })


@router.delete("/cancel/{job_id}")
async def cancel_simulation(job_id: str):
    """Cancel a running or pending simulation."""
    if job_id in ACTIVE_SIMULATIONS:
        job = ACTIVE_SIMULATIONS[job_id]
        if job["status"] in ["pending", "queued"]:
            job["status"] = "cancelled"
            job["completed_at"] = datetime.utcnow().isoformat() + "Z"
            save_job(job_id, job)
            return JSONResponse({
                "success": True,
                "job_id": job_id,
                "message": "Simulation cancelled",
            })
        elif job["status"] == "running":
            return JSONResponse({
                "success": False,
                "job_id": job_id,
                "message": "Cannot cancel running simulation",
            }, status_code=400)
        else:
            return JSONResponse({
                "success": False,
                "job_id": job_id,
                "message": f"Simulation already {job['status']}",
            }, status_code=400)

    raise HTTPException(status_code=404, detail=f"Job {job_id} not found")


@router.get("/presets")
async def list_simulation_presets():
    """List predefined simulation presets matching CLI examples."""
    presets = [
        {
            "code": "example1_monolayer",
            "name": "Example 1: Monolayer LDPE",
            "description": "Simple monolayer LDPE with BHT antioxidant",
            "source": "examples/example1.py",
            "config": {
                "layers": [
                    {"index": 1, "polymer": "LDPE", "thickness": 100, "thickness_unit": "um", "C0": 1000}
                ],
                "steps": [
                    {"index": 1, "temperature_C": 40, "duration": 10, "duration_unit": "days"}
                ],
                "food": {"category": "fatty"},
            },
        },
        {
            "code": "example2_bilayer",
            "name": "Example 2: Bilayer PET/LDPE",
            "description": "Functional barrier with glassy PET layer",
            "source": "examples/example2.py",
            "config": {
                "layers": [
                    {"index": 1, "polymer": "gPET", "thickness": 20, "thickness_unit": "um", "C0": 0},
                    {"index": 2, "polymer": "LDPE", "thickness": 100, "thickness_unit": "um", "C0": 1000},
                ],
                "steps": [
                    {"index": 1, "temperature_C": 40, "duration": 10, "duration_unit": "days"}
                ],
                "food": {"category": "fatty"},
            },
        },
        {
            "code": "example3_trilayer",
            "name": "Example 3: Trilayer ABA",
            "description": "Symmetric ABA structure for set-off studies",
            "source": "examples/example3.py",
            "config": {
                "layers": [
                    {"index": 1, "polymer": "LDPE", "thickness": 50, "thickness_unit": "um", "C0": 0},
                    {"index": 2, "polymer": "PP", "thickness": 200, "thickness_unit": "um", "C0": 1000},
                    {"index": 3, "polymer": "LDPE", "thickness": 50, "thickness_unit": "um", "C0": 0},
                ],
                "steps": [
                    {"index": 1, "temperature_C": 25, "duration": 90, "duration_unit": "days", "with_food": False},
                    {"index": 2, "temperature_C": 40, "duration": 10, "duration_unit": "days", "with_food": True},
                ],
                "food": {"category": "fatty"},
            },
        },
        {
            "code": "hotfill_scenario",
            "name": "Hot Fill Scenario",
            "description": "Hot filling followed by storage",
            "source": "custom",
            "config": {
                "layers": [
                    {"index": 1, "polymer": "PET", "thickness": 300, "thickness_unit": "um", "C0": 500}
                ],
                "steps": [
                    {"index": 1, "temperature_C": 85, "duration": 30, "duration_unit": "min"},
                    {"index": 2, "temperature_C": 25, "duration": 6, "duration_unit": "months"},
                ],
                "food": {"category": "aqueous"},
            },
        },
    ]

    return JSONResponse({
        "success": True,
        "presets": presets,
        "count": len(presets),
    })


@router.post("/validate")
async def validate_simulation_config(config: SimulationInput):
    """Validate simulation configuration before running."""
    errors = []
    warnings = []

    # Validate layers
    if not config.layers:
        errors.append("At least one layer is required")
    else:
        indices = [l.index for l in config.layers]
        if sorted(indices) != list(range(1, len(indices) + 1)):
            errors.append("Layer indices must be consecutive starting from 1")

        # Check for missing D values
        for layer in config.layers:
            if layer.D is None and (layer.D_override is None or layer.D_override.use_computed):
                warnings.append(f"Layer {layer.index}: D will be estimated using Piringer model")

    # Validate steps
    if not config.steps:
        errors.append("At least one step is required")
    else:
        step_indices = [s.index for s in config.steps]
        if sorted(step_indices) != list(range(1, len(step_indices) + 1)):
            errors.append("Step indices must be consecutive starting from 1")

    # Calculate totals
    total_duration_s = sum(
        convert_duration(s.duration, s.duration_unit)
        for s in config.steps
    )

    return JSONResponse({
        "success": len(errors) == 0,
        "valid": len(errors) == 0,
        "errors": errors,
        "warnings": warnings,
        "layer_count": len(config.layers),
        "step_count": len(config.steps),
        "total_duration_days": total_duration_s / 86400,
    })
