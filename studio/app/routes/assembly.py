"""
Assembly Routes - Material Layer Management

API endpoints for defining multilayer packaging structures.
"""

import sys
from pathlib import Path
from typing import List, Optional

from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

import yaml

# Add parent paths
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

router = APIRouter()

# ========== DATA MODELS ==========

class ParameterOverride(BaseModel):
    use_computed: bool = True
    override_value: Optional[float] = None


class LayerInput(BaseModel):
    index: int = Field(..., ge=1, le=10, description="Layer index (1 = food contact)")
    polymer: str = Field(..., description="Polymer code (e.g., 'LDPE', 'PET')")
    thickness: float = Field(..., gt=0, description="Thickness value")
    thickness_unit: str = Field(default="um", description="Thickness unit")
    temperature_C: float = Field(default=25.0, description="Temperature in Celsius")
    is_glassy: bool = Field(default=False, description="Is polymer glassy (functional barrier)")
    is_plasticized: bool = Field(default=False, description="Is polymer plasticized")
    D_override: Optional[ParameterOverride] = None
    k_override: Optional[ParameterOverride] = None
    k0_override: Optional[ParameterOverride] = None


class AssemblyInput(BaseModel):
    name: str = Field(default="Untitled", description="Assembly name")
    layers: List[LayerInput] = Field(..., min_length=1, max_length=10)


class PolymerInfo(BaseModel):
    code: str
    name: str
    Tg: float  # Glass transition temperature
    density: float  # kg/m³
    icon: str
    description: Optional[str] = None


# ========== POLYMER DATABASE ==========

def load_polymers() -> dict:
    """Load polymer database from YAML."""
    data_dir = Path(__file__).parent.parent.parent / "data"
    polymers_file = data_dir / "polymers.yaml"

    if polymers_file.exists():
        with open(polymers_file, 'r') as f:
            data = yaml.safe_load(f)
            return {p['code']: p for p in data.get('polymers', [])}

    # Default polymers if file doesn't exist
    return {
        "LDPE": {"code": "LDPE", "name": "Low Density Polyethylene", "Tg": -120, "density": 920, "icon": "♳"},
        "HDPE": {"code": "HDPE", "name": "High Density Polyethylene", "Tg": -120, "density": 960, "icon": "♴"},
        "PP": {"code": "PP", "name": "Polypropylene", "Tg": -10, "density": 905, "icon": "♵"},
        "PET": {"code": "PET", "name": "Polyethylene Terephthalate", "Tg": 75, "density": 1380, "icon": "♳"},
        "gPET": {"code": "gPET", "name": "Glassy PET", "Tg": 75, "density": 1380, "icon": "♳", "is_glassy": True},
        "wPET": {"code": "wPET", "name": "Wet/Plasticized PET", "Tg": 40, "density": 1380, "icon": "♳"},
        "PS": {"code": "PS", "name": "Polystyrene", "Tg": 100, "density": 1050, "icon": "♶"},
        "PA6": {"code": "PA6", "name": "Polyamide 6 (Nylon)", "Tg": 50, "density": 1130, "icon": "♷"},
        "PVDC": {"code": "PVDC", "name": "Polyvinylidene Chloride", "Tg": -17, "density": 1700, "icon": "♸"},
        "EVOH": {"code": "EVOH", "name": "Ethylene Vinyl Alcohol", "Tg": 55, "density": 1140, "icon": "♹"},
    }


POLYMERS_DB = load_polymers()


# ========== MATERIALS DISCOVERY ==========

def discover_materials_from_layer() -> dict:
    """
    Discover all available materials from patankar.layer module.

    Returns dict of material_code -> material_info including:
    - code, name, category, Tg, Tm, default_thickness, etc.
    """
    try:
        from patankar import layer as layer_module
        import inspect

        # Get the base layer class
        base_layer = getattr(layer_module, 'layer', None)
        if base_layer is None:
            return {}

        materials = {}

        # Iterate through all classes in layer module
        for name, cls in inspect.getmembers(layer_module, inspect.isclass):
            # Skip the base class and non-layer classes
            if name == 'layer' or not issubclass(cls, base_layer):
                continue

            try:
                # Create a temporary instance to get properties
                instance = cls()

                # Determine category based on class name and properties
                category = "polymer"
                if "Adhesive" in name:
                    category = "adhesive"
                elif name in ["Paper", "Cardboard"]:
                    category = "paper"
                elif name == "air":
                    category = "gas"

                # Get Tg safely
                tg = None
                try:
                    tg_prop = instance.Tg
                    if isinstance(tg_prop, tuple):
                        tg = tg_prop[0]
                    elif isinstance(tg_prop, (int, float)):
                        tg = tg_prop
                except:
                    pass

                # Get Tm safely
                tm = None
                try:
                    tm_prop = instance.Tm
                    if tm_prop is not None:
                        if isinstance(tm_prop, tuple):
                            tm = tm_prop[0]
                        elif isinstance(tm_prop, (int, float)):
                            tm = tm_prop
                except:
                    pass

                # Get default thickness in micrometers
                default_l = 100  # default
                try:
                    l_val = instance.l[0] if hasattr(instance.l, '__iter__') else instance.l
                    default_l = l_val * 1e6  # Convert from m to µm
                except:
                    pass

                # Get material description
                material_name = getattr(instance, 'layermaterial', name)
                layer_code = getattr(instance, 'layercode', name)

                # Determine if it's a barrier (high Tg)
                is_barrier = tg is not None and tg > 50

                materials[name] = {
                    "code": name,
                    "layer_code": layer_code,
                    "name": material_name if isinstance(material_name, str) else name,
                    "category": category,
                    "Tg": tg,
                    "Tm": tm,
                    "default_thickness_um": round(default_l, 1),
                    "is_barrier": is_barrier,
                    "description": cls.__doc__.strip().split('\n')[0] if cls.__doc__ else None,
                }

            except Exception as e:
                # Skip materials that can't be instantiated
                continue

        return materials

    except ImportError:
        return {}


# Cache discovered materials
_DISCOVERED_MATERIALS = None


def get_all_materials() -> dict:
    """Get all materials from both YAML and layer.py discovery."""
    global _DISCOVERED_MATERIALS

    if _DISCOVERED_MATERIALS is None:
        _DISCOVERED_MATERIALS = discover_materials_from_layer()

    # Merge with POLYMERS_DB (YAML takes precedence for icons/descriptions)
    all_materials = dict(_DISCOVERED_MATERIALS)
    for code, info in POLYMERS_DB.items():
        if code in all_materials:
            all_materials[code].update(info)
        else:
            all_materials[code] = info

    return all_materials


# ========== ESTIMATION FUNCTIONS ==========

def estimate_diffusivity(polymer: str, mw: float, temp_C: float) -> float:
    """
    Estimate diffusivity using Piringer model.

    D = D0 * exp(-Ap/T) * exp(-B * MW^n)

    Returns D in m²/s
    """
    # Piringer model parameters (simplified)
    Ap_values = {
        "LDPE": 0, "HDPE": 0, "PP": 0,
        "PET": 1577, "gPET": 1577, "wPET": 1000,
        "PS": 0, "PA6": 0, "PVDC": 0, "EVOH": 0,
    }

    Ap = Ap_values.get(polymer, 0)
    T_K = temp_C + 273.15

    # Simplified Piringer equation
    # log10(D) = Ap - 10454/T - 0.003 * MW
    if Ap > 0:
        log_D = Ap - 10454 / T_K - 0.003 * mw
    else:
        # For polyolefins: log10(D) = -4 - 10454/T - 0.003 * MW
        log_D = -4 - 10454 / T_K - 0.003 * mw

    D = 10 ** log_D * 1e-4  # Convert from cm²/s to m²/s
    return max(D, 1e-20)  # Ensure positive


def estimate_partition(polymer: str, logP: float, temp_C: float) -> float:
    """
    Estimate partition coefficient using FHP model.

    Returns k (dimensionless)
    """
    # Simplified: k ~ 1 for most cases, adjusted by logP
    # High logP (lipophilic) → higher k in lipophilic polymers
    if polymer in ["LDPE", "HDPE", "PP"]:
        k = 10 ** (0.1 * logP - 0.5) if logP > 5 else 1.0
    else:
        k = 1.0

    return max(k, 0.01)


# ========== API ENDPOINTS ==========

@router.get("/polymers")
async def list_polymers():
    """List all available polymers."""
    return JSONResponse({
        "success": True,
        "polymers": list(POLYMERS_DB.values()),
        "count": len(POLYMERS_DB),
    })


@router.get("/materials")
async def list_all_materials(category: str = None):
    """
    List all available materials discovered from layer.py.

    Categories: polymer, adhesive, paper, gas

    This endpoint dynamically discovers all material classes from the
    patankar.layer module, including their Tg, Tm, and default properties.
    """
    all_materials = get_all_materials()

    # Filter by category if specified
    if category:
        all_materials = {k: v for k, v in all_materials.items()
                         if v.get("category") == category}

    # Sort by category and then by name
    sorted_materials = sorted(
        all_materials.values(),
        key=lambda x: (x.get("category", "z"), x.get("name", ""))
    )

    return JSONResponse({
        "success": True,
        "materials": sorted_materials,
        "count": len(sorted_materials),
        "categories": ["polymer", "adhesive", "paper", "gas"],
    })


@router.get("/materials/{code}")
async def get_material(code: str):
    """Get details for a specific material."""
    all_materials = get_all_materials()

    if code not in all_materials:
        raise HTTPException(status_code=404, detail=f"Material '{code}' not found")

    return JSONResponse({
        "success": True,
        "material": all_materials[code],
    })


@router.get("/materials/{code}/tg")
async def get_material_tg(code: str):
    """Get glass transition temperature for a material."""
    all_materials = get_all_materials()

    if code not in all_materials:
        raise HTTPException(status_code=404, detail=f"Material '{code}' not found")

    material = all_materials[code]
    tg = material.get("Tg")

    return JSONResponse({
        "success": True,
        "code": code,
        "Tg": tg,
        "Tg_unit": "°C",
        "is_barrier": material.get("is_barrier", False),
    })


@router.get("/materials/{code}/density")
async def get_material_density(code: str, temperature_C: float = 25.0):
    """
    Get density for a material at specified temperature.

    Uses patankar layer.density(T) method which incorporates temperature dependence:
        e.g., LDPE.density(T) = 920 * (1 - 3*(T - Td) * 20e-5)

    Fallback chain:
    1. patankar layer.density(T) method (T-dependent)
    2. patankar layer.rho attribute (may not be T-dependent)
    3. YAML/default data

    Args:
        code: Material code (e.g., 'LDPE', 'PP', 'PET')
        temperature_C: Temperature in Celsius (default 25°C)

    Returns:
        Density in kg/m³ and related metadata
    """
    all_materials = get_all_materials()

    if code not in all_materials:
        raise HTTPException(status_code=404, detail=f"Material '{code}' not found")

    # Get density from patankar layer instance (with temperature support)
    rho = None
    rho_unit = "kg/m³"
    source = "yaml"  # Default source
    is_temperature_dependent = False

    try:
        from patankar import layer as layer_module

        # Get the layer class
        LayerClass = getattr(layer_module, code, None)
        if LayerClass is not None:
            # Instantiate with temperature
            instance = LayerClass(l=(100, "um"), T=(temperature_C, "degC"))

            # Priority 1: Use density(T) method which has T-dependence
            if hasattr(instance, 'density') and callable(instance.density):
                try:
                    rho_result = instance.density(temperature_C)
                    # density() returns tuple (array([920.]), "kg/m**3")
                    if isinstance(rho_result, tuple):
                        rho_val = rho_result[0]
                        # Handle numpy array: array([920.]) -> 920.0
                        if hasattr(rho_val, '__iter__') and not isinstance(rho_val, str):
                            rho = float(rho_val[0])
                        else:
                            rho = float(rho_val)
                        rho_unit = rho_result[1] if len(rho_result) > 1 else "kg/m³"
                    else:
                        rho = float(rho_result)
                    source = "patankar"
                    is_temperature_dependent = True
                except Exception:
                    pass

            # Priority 2: Fallback to rho attribute (may not be T-dependent)
            if rho is None and hasattr(instance, 'rho') and instance.rho is not None:
                rho_val = instance.rho
                # Handle numpy array
                if hasattr(rho_val, '__iter__') and not isinstance(rho_val, str):
                    rho = float(rho_val[0])
                else:
                    rho = float(rho_val)
                source = "patankar"
                is_temperature_dependent = False

            if hasattr(instance, 'rhounit'):
                rho_unit = instance.rhounit
    except Exception as e:
        # Fall back to YAML data if patankar fails
        pass

    # Priority 3: Fallback to YAML/default data
    if rho is None:
        material = all_materials[code]
        rho = material.get("density")
        source = "yaml"
        is_temperature_dependent = False

    # Priority 4: Hardcoded defaults for common materials
    if rho is None:
        DEFAULT_DENSITIES = {
            'LDPE': 920, 'LLDPE': 920, 'HDPE': 960,
            'PP': 910, 'gPET': 1380, 'wPET': 1380, 'PET': 1380, 'rPET': 1380,
            'PS': 1050, 'PMMA': 1190, 'PA6': 1130, 'PA66': 1140,
            'PBT': 1310, 'PEN': 1360, 'HIPS': 1040, 'SBS': 940,
        }
        rho = DEFAULT_DENSITIES.get(code)
        source = "default"
        is_temperature_dependent = False

    if rho is None:
        raise HTTPException(status_code=404, detail=f"Density not available for material '{code}'")

    return JSONResponse({
        "success": True,
        "code": code,
        "rho": rho,
        "rho_unit": rho_unit,
        "temperature_C": temperature_C,
        "source": source,  # 'patankar', 'yaml', or 'default'
        "is_temperature_dependent": is_temperature_dependent,
    })


@router.get("/polymers/{code}")
async def get_polymer(code: str):
    """Get details for a specific polymer."""
    if code not in POLYMERS_DB:
        raise HTTPException(status_code=404, detail=f"Polymer '{code}' not found")

    return JSONResponse({
        "success": True,
        "polymer": POLYMERS_DB[code],
    })


@router.post("/validate")
async def validate_assembly(assembly: AssemblyInput):
    """
    Validate an assembly configuration.

    Checks:
    - All polymers exist
    - Layer indices are correct
    - Thicknesses are reasonable
    """
    errors = []
    warnings = []

    # Check layer indices
    indices = [layer.index for layer in assembly.layers]
    if sorted(indices) != list(range(1, len(indices) + 1)):
        errors.append("Layer indices must be consecutive starting from 1")

    # Check polymers
    for layer in assembly.layers:
        if layer.polymer not in POLYMERS_DB:
            errors.append(f"Unknown polymer: {layer.polymer}")

        # Check thickness (convert to meters for validation)
        thickness_m = convert_thickness(layer.thickness, layer.thickness_unit)
        if thickness_m < 1e-9:
            warnings.append(f"Layer {layer.index}: Very thin layer ({layer.thickness} {layer.thickness_unit})")
        if thickness_m > 0.1:
            warnings.append(f"Layer {layer.index}: Very thick layer ({layer.thickness} {layer.thickness_unit})")

    # Calculate total thickness
    total_thickness_m = sum(
        convert_thickness(layer.thickness, layer.thickness_unit)
        for layer in assembly.layers
    )

    return JSONResponse({
        "success": len(errors) == 0,
        "valid": len(errors) == 0,
        "errors": errors,
        "warnings": warnings,
        "total_thickness_m": total_thickness_m,
        "total_thickness_um": total_thickness_m * 1e6,
    })


@router.post("/estimate")
async def estimate_properties(assembly: AssemblyInput, mw: float = 500.0, logP: float = 5.0):
    """
    Estimate D and k for each layer based on substance properties.

    Args:
        assembly: The layer assembly
        mw: Molecular weight of substance (g/mol)
        logP: Log partition coefficient (octanol/water)
    """
    results = []

    for layer in assembly.layers:
        D_computed = estimate_diffusivity(layer.polymer, mw, layer.temperature_C)
        k_computed = estimate_partition(layer.polymer, logP, layer.temperature_C)

        # Apply overrides if specified
        D_final = D_computed
        if layer.D_override and not layer.D_override.use_computed:
            D_final = layer.D_override.override_value or D_computed

        k_final = k_computed
        if layer.k_override and not layer.k_override.use_computed:
            k_final = layer.k_override.override_value or k_computed

        results.append({
            "layer_index": layer.index,
            "polymer": layer.polymer,
            "temperature_C": layer.temperature_C,
            "D_computed": D_computed,
            "D_final": D_final,
            "D_override_applied": layer.D_override is not None and not layer.D_override.use_computed,
            "k_computed": k_computed,
            "k_final": k_final,
            "k_override_applied": layer.k_override is not None and not layer.k_override.use_computed,
        })

    return JSONResponse({
        "success": True,
        "mw": mw,
        "logP": logP,
        "layers": results,
    })


def convert_thickness(value: float, unit: str) -> float:
    """Convert thickness to meters."""
    conversions = {
        "nm": 1e-9,
        "um": 1e-6,
        "µm": 1e-6,
        "mm": 1e-3,
        "cm": 1e-2,
        "m": 1.0,
        "in": 0.0254,
    }
    return value * conversions.get(unit.lower(), 1e-6)
