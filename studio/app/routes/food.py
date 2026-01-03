"""
Food & Conditions Routes

API endpoints for food types, geometry, and contact conditions.
"""

import sys
import math
from pathlib import Path
from typing import List, Optional

from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

import yaml

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

router = APIRouter()


# ========== DATA MODELS ==========

class GeometryInput(BaseModel):
    shape: str = Field(..., description="Shape type: cylinder, bottle, rectangle, sphere")
    dimensions: dict = Field(..., description="Shape-specific dimensions in mm")


class ContactStepInput(BaseModel):
    index: int = Field(..., ge=1, le=10)
    condition_type: str = Field(default="storage", description="Type: storage, hotfill, transport")
    temperature_C: float = Field(default=25.0)
    duration: float = Field(..., gt=0)
    duration_unit: str = Field(default="days")
    with_food: bool = Field(default=True, description="False for set-off (no food contact)")


class FoodInput(BaseModel):
    category: str = Field(..., description="Food category: fatty, aqueous, acidic, alcoholic, dry")
    texture: str = Field(default="liquid", description="Texture: liquid, semisolid, solid")
    process: str = Field(default="ambient", description="Process: ambient, hotfilled, frozen")
    simulant: Optional[str] = Field(None, description="Food simulant code")
    use_real_food: bool = Field(default=False)


class ScenarioInput(BaseModel):
    food: FoodInput
    geometry: GeometryInput
    steps: List[ContactStepInput] = Field(..., min_length=1, max_length=10)


# ========== DATABASE LOADERS ==========

def load_foods() -> dict:
    """Load food types database."""
    return {
        "fatty": {
            "code": "fatty",
            "name": "Fatty foods",
            "description": "Foods with fat content > 20%",
            "default_simulant": "olive_oil",
            "icon": "🥩",
            "examples": ["butter", "cheese", "meat", "nuts"],
        },
        "aqueous": {
            "code": "aqueous",
            "name": "Aqueous foods",
            "description": "Water-based foods, pH > 4.5",
            "default_simulant": "water",
            "icon": "💧",
            "examples": ["milk", "juices", "soft drinks"],
        },
        "acidic": {
            "code": "acidic",
            "name": "Acidic foods",
            "description": "Foods with pH < 4.5",
            "default_simulant": "acetic_acid_3",
            "icon": "🍋",
            "examples": ["vinegar", "citrus juices", "pickles"],
        },
        "alcoholic": {
            "code": "alcoholic",
            "name": "Alcoholic foods",
            "description": "Foods containing alcohol > 5%",
            "default_simulant": "ethanol_20",
            "icon": "🍷",
            "examples": ["wine", "beer", "spirits"],
        },
        "dry": {
            "code": "dry",
            "name": "Dry foods",
            "description": "Foods with water activity < 0.7",
            "default_simulant": "tenax",
            "icon": "🍞",
            "examples": ["bread", "cereals", "pasta"],
        },
    }


def load_simulants() -> dict:
    """Load food simulants database (EU 10/2011) with density data."""
    return {
        "water": {
            "code": "water",
            "name": "Simulant A - Water",
            "description": "For aqueous foods",
            "composition": "Distilled water",
            "category": "aqueous",
            "density_kg_m3": 997,  # at 25°C
            "density_temp_coeff": -0.0003,  # dρ/dT per °C
        },
        "acetic_acid_3": {
            "code": "acetic_acid_3",
            "name": "Simulant B - Acetic acid 3%",
            "description": "For acidic foods",
            "composition": "3% (w/v) acetic acid in water",
            "category": "acidic",
            "density_kg_m3": 1006,  # at 25°C
            "density_temp_coeff": -0.0003,
        },
        "ethanol_10": {
            "code": "ethanol_10",
            "name": "Simulant C - Ethanol 10%",
            "description": "For alcoholic foods < 10%",
            "composition": "10% (v/v) ethanol in water",
            "category": "alcoholic",
            "density_kg_m3": 982,  # at 25°C
            "density_temp_coeff": -0.0008,
        },
        "ethanol_20": {
            "code": "ethanol_20",
            "name": "Simulant C - Ethanol 20%",
            "description": "For alcoholic foods 10-20%",
            "composition": "20% (v/v) ethanol in water",
            "category": "alcoholic",
            "density_kg_m3": 969,  # at 25°C
            "density_temp_coeff": -0.0009,
        },
        "ethanol_50": {
            "code": "ethanol_50",
            "name": "Simulant D1 - Ethanol 50%",
            "description": "For fatty foods (alternative to oil)",
            "composition": "50% (v/v) ethanol in water",
            "category": "fatty",
            "density_kg_m3": 914,  # at 25°C
            "density_temp_coeff": -0.0010,
        },
        "olive_oil": {
            "code": "olive_oil",
            "name": "Simulant D2 - Olive oil",
            "description": "For fatty foods (reference)",
            "composition": "Rectified olive oil",
            "category": "fatty",
            "density_kg_m3": 920,  # at 25°C
            "density_temp_coeff": -0.0007,
        },
        "tenax": {
            "code": "tenax",
            "name": "Simulant E - Tenax",
            "description": "For dry foods",
            "composition": "Poly(2,6-diphenyl-p-phenylene oxide)",
            "category": "dry",
            "density_kg_m3": None,  # Solid adsorbent - no liquid density
            "density_temp_coeff": None,
        },
        "ethanol_95": {
            "code": "ethanol_95",
            "name": "Ethanol 95%",
            "description": "For high-fat and worst-case testing",
            "composition": "95% (v/v) ethanol in water",
            "category": "fatty",
            "density_kg_m3": 789,  # at 25°C
            "density_temp_coeff": -0.0011,
        },
        # Additional simulants for patankar compatibility
        "ethanol": {
            "code": "ethanol",
            "name": "Ethanol (generic)",
            "description": "Generic ethanol simulant",
            "composition": "Ethanol solution",
            "category": "fatty",
            "density_kg_m3": 914,  # Default to 50% ethanol density
            "density_temp_coeff": -0.0010,
        },
        "isooctane": {
            "code": "isooctane",
            "name": "Isooctane",
            "description": "Alternative fatty simulant",
            "composition": "2,2,4-trimethylpentane",
            "category": "fatty",
            "density_kg_m3": 692,  # at 25°C
            "density_temp_coeff": -0.0012,
        },
    }


FOODS_DB = load_foods()
SIMULANTS_DB = load_simulants()


# ========== GEOMETRY CALCULATIONS ==========

def calculate_geometry(shape: str, dimensions: dict) -> dict:
    """
    Calculate volume and surface area from shape and dimensions.

    All dimensions expected in mm, returns m³ and m².
    """
    if shape == "cylinder":
        # dimensions: radius, height (mm)
        r = dimensions.get("radius", 50) / 1000  # Convert to m
        h = dimensions.get("height", 100) / 1000
        volume = math.pi * r**2 * h
        surface = 2 * math.pi * r * h + 2 * math.pi * r**2  # Lateral + top + bottom
        # For food contact, often only inner surface matters
        inner_surface = 2 * math.pi * r * h + math.pi * r**2  # Lateral + bottom

    elif shape == "bottle":
        # dimensions: body_radius, body_height, neck_radius, neck_height (mm)
        rb = dimensions.get("body_radius", 40) / 1000
        hb = dimensions.get("body_height", 200) / 1000
        rn = dimensions.get("neck_radius", 15) / 1000
        hn = dimensions.get("neck_height", 50) / 1000
        volume = math.pi * rb**2 * hb + math.pi * rn**2 * hn
        inner_surface = 2 * math.pi * rb * hb + math.pi * rb**2 + 2 * math.pi * rn * hn
        surface = inner_surface

    elif shape == "rectangle" or shape == "tray":
        # dimensions: length, width, height (mm)
        l = dimensions.get("length", 200) / 1000
        w = dimensions.get("width", 150) / 1000
        h = dimensions.get("height", 50) / 1000
        volume = l * w * h
        surface = 2 * (l * w + l * h + w * h)
        inner_surface = l * w + 2 * l * h + 2 * w * h  # Bottom + sides

    elif shape == "sphere":
        # dimensions: radius (mm)
        r = dimensions.get("radius", 50) / 1000
        volume = (4/3) * math.pi * r**3
        surface = 4 * math.pi * r**2
        inner_surface = surface

    elif shape == "pouch":
        # dimensions: width, height (mm) - flat pouch
        w = dimensions.get("width", 150) / 1000
        h = dimensions.get("height", 200) / 1000
        thickness = dimensions.get("thickness", 10) / 1000  # Pouch thickness when filled
        volume = w * h * thickness
        surface = 2 * w * h  # Both sides
        inner_surface = surface

    else:
        raise ValueError(f"Unknown shape: {shape}")

    return {
        "volume_m3": volume,
        "surface_m2": surface,
        "inner_surface_m2": inner_surface,
        "volume_cm3": volume * 1e6,
        "surface_cm2": surface * 1e4,
        "vs_ratio_cm": (volume * 1e6) / (inner_surface * 1e4) if inner_surface > 0 else 0,
    }


def convert_duration(value: float, unit: str) -> float:
    """Convert duration to seconds."""
    conversions = {
        "s": 1,
        "sec": 1,
        "seconds": 1,
        "min": 60,
        "minutes": 60,
        "h": 3600,
        "hours": 3600,
        "days": 86400,
        "d": 86400,
        "weeks": 604800,
        "w": 604800,
        "months": 2592000,  # 30 days
        "years": 31536000,  # 365 days
    }
    return value * conversions.get(unit.lower(), 86400)


# ========== API ENDPOINTS ==========

@router.get("/categories")
async def list_food_categories():
    """List all food categories."""
    return JSONResponse({
        "success": True,
        "categories": list(FOODS_DB.values()),
        "count": len(FOODS_DB),
    })


@router.get("/simulants")
async def list_simulants():
    """List all food simulants."""
    return JSONResponse({
        "success": True,
        "simulants": list(SIMULANTS_DB.values()),
        "count": len(SIMULANTS_DB),
    })


@router.get("/simulants/{category}")
async def get_simulants_for_category(category: str):
    """Get recommended simulants for a food category."""
    if category not in FOODS_DB:
        raise HTTPException(status_code=404, detail=f"Unknown category: {category}")

    food = FOODS_DB[category]
    simulants = [s for s in SIMULANTS_DB.values() if s.get("category") == category]

    return JSONResponse({
        "success": True,
        "category": category,
        "default_simulant": food.get("default_simulant"),
        "simulants": simulants,
    })


@router.get("/simulants/density/{code}")
async def get_simulant_density(code: str, temperature_C: float = 25.0):
    """
    Get density for a food simulant at specified temperature.

    Supports temperature-dependent density using linear model:
    ρ(T) = ρ_25 * (1 + α * (T - 25))

    Args:
        code: Simulant code (e.g., 'water', 'ethanol_50', 'olive_oil')
        temperature_C: Temperature in Celsius (default 25°C)

    Returns:
        Density in kg/m³ and related metadata
    """
    if code not in SIMULANTS_DB:
        raise HTTPException(status_code=404, detail=f"Unknown simulant: {code}")

    simulant = SIMULANTS_DB[code]
    rho_25 = simulant.get("density_kg_m3")

    if rho_25 is None:
        raise HTTPException(
            status_code=400,
            detail=f"Density not available for simulant '{code}' (solid adsorbent)"
        )

    # Apply temperature correction
    temp_coeff = simulant.get("density_temp_coeff", 0) or 0
    rho_T = rho_25 * (1 + temp_coeff * (temperature_C - 25.0))

    return JSONResponse({
        "success": True,
        "code": code,
        "density_kg_m3": round(rho_T, 2),
        "density_25C_kg_m3": rho_25,
        "temperature_C": temperature_C,
        "temp_coeff": temp_coeff,
        "unit": "kg/m³",
    })


@router.get("/food/density")
async def get_food_density(
    category: str = "fatty",
    simulant: Optional[str] = None,
    temperature_C: float = 25.0
):
    """
    Get density for a food category or specific simulant.

    If simulant is provided, returns that simulant's density.
    Otherwise, returns the default simulant density for the category.

    Args:
        category: Food category (fatty, aqueous, etc.)
        simulant: Specific simulant code (optional)
        temperature_C: Temperature in Celsius (default 25°C)

    Returns:
        Density in kg/m³ with source information
    """
    # Try to get density from patankar food instance first
    rho = None
    source = "database"

    # Determine which simulant to use
    if simulant and simulant in SIMULANTS_DB:
        sim = SIMULANTS_DB[simulant]
    elif category in FOODS_DB:
        default_sim = FOODS_DB[category].get("default_simulant")
        sim = SIMULANTS_DB.get(default_sim) if default_sim else None
    else:
        sim = None

    if sim:
        rho_25 = sim.get("density_kg_m3")
        if rho_25 is not None:
            temp_coeff = sim.get("density_temp_coeff", 0) or 0
            rho = rho_25 * (1 + temp_coeff * (temperature_C - 25.0))
            source = f"simulant:{sim['code']}"

    # Try patankar food module as fallback
    if rho is None:
        try:
            from patankar import food as food_module

            # Create a food instance to get density
            class TestFood(food_module.realfood, food_module.liquid, food_module.fat):
                pass

            food_instance = TestFood(contacttemperature=(temperature_C, "degC"))
            if hasattr(food_instance, 'density') and food_instance.density is not None:
                rho = float(food_instance.density)
                source = "patankar"
        except Exception:
            pass

    # Default fallback
    if rho is None:
        rho = 1000.0  # Default water density
        source = "default"

    return JSONResponse({
        "success": True,
        "category": category,
        "simulant": simulant,
        "density_kg_m3": round(rho, 2),
        "temperature_C": temperature_C,
        "source": source,
        "unit": "kg/m³",
    })


@router.post("/geometry/calculate")
async def calculate_geometry_endpoint(geometry: GeometryInput):
    """Calculate volume and surface area from geometry."""
    try:
        result = calculate_geometry(geometry.shape, geometry.dimensions)
        return JSONResponse({
            "success": True,
            "shape": geometry.shape,
            "dimensions": geometry.dimensions,
            **result,
        })
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.get("/shapes")
async def list_shapes():
    """List available packaging shapes."""
    shapes = [
        {
            "code": "cylinder",
            "name": "Cylinder",
            "icon": "🥫",
            "dimensions": ["radius", "height"],
            "description": "Cylindrical container (cans, jars)",
        },
        {
            "code": "bottle",
            "name": "Bottle",
            "icon": "🍼",
            "dimensions": ["body_radius", "body_height", "neck_radius", "neck_height"],
            "description": "Bottle with neck",
        },
        {
            "code": "rectangle",
            "name": "Rectangle/Tray",
            "icon": "🍱",
            "dimensions": ["length", "width", "height"],
            "description": "Rectangular container or tray",
        },
        {
            "code": "pouch",
            "name": "Pouch/Sachet",
            "icon": "🛍️",
            "dimensions": ["width", "height", "thickness"],
            "description": "Flexible pouch",
        },
        {
            "code": "sphere",
            "name": "Sphere",
            "icon": "⚽",
            "dimensions": ["radius"],
            "description": "Spherical container",
        },
    ]
    return JSONResponse({
        "success": True,
        "shapes": shapes,
    })


@router.post("/scenario/validate")
async def validate_scenario(scenario: ScenarioInput):
    """Validate a complete scenario configuration."""
    errors = []
    warnings = []

    # Validate food
    if scenario.food.category not in FOODS_DB:
        errors.append(f"Unknown food category: {scenario.food.category}")

    if scenario.food.simulant and scenario.food.simulant not in SIMULANTS_DB:
        errors.append(f"Unknown simulant: {scenario.food.simulant}")

    # Validate geometry
    try:
        geom = calculate_geometry(scenario.geometry.shape, scenario.geometry.dimensions)
        if geom["volume_m3"] < 1e-9:
            warnings.append("Very small volume - check dimensions")
        if geom["vs_ratio_cm"] > 10:
            warnings.append("High V/S ratio - migration may be slow")
    except ValueError as e:
        errors.append(str(e))
        geom = None

    # Validate steps
    step_indices = [s.index for s in scenario.steps]
    if sorted(step_indices) != list(range(1, len(step_indices) + 1)):
        errors.append("Step indices must be consecutive starting from 1")

    total_duration_s = sum(
        convert_duration(s.duration, s.duration_unit)
        for s in scenario.steps
    )

    # Check for set-off before food contact
    has_setoff = any(not s.with_food for s in scenario.steps)
    if has_setoff:
        first_food_step = next((s for s in scenario.steps if s.with_food), None)
        if first_food_step and first_food_step.index == 1:
            warnings.append("First step has food contact - set-off steps should come first")

    return JSONResponse({
        "success": len(errors) == 0,
        "valid": len(errors) == 0,
        "errors": errors,
        "warnings": warnings,
        "geometry": geom,
        "total_duration_s": total_duration_s,
        "total_duration_days": total_duration_s / 86400,
        "has_setoff": has_setoff,
        "step_count": len(scenario.steps),
    })


@router.get("/conditions/presets")
async def list_condition_presets():
    """List common condition presets."""
    presets = [
        {
            "code": "eu_standard",
            "name": "EU Standard Test",
            "description": "10 days at 40°C (EU 10/2011)",
            "steps": [{"temperature_C": 40, "duration": 10, "unit": "days"}],
        },
        {
            "code": "accelerated",
            "name": "Accelerated Test",
            "description": "10 days at 60°C",
            "steps": [{"temperature_C": 60, "duration": 10, "unit": "days"}],
        },
        {
            "code": "refrigerated",
            "name": "Refrigerated Storage",
            "description": "Long-term at 4°C",
            "steps": [{"temperature_C": 4, "duration": 6, "unit": "months"}],
        },
        {
            "code": "hotfill",
            "name": "Hot Fill Process",
            "description": "Hot fill + storage",
            "steps": [
                {"temperature_C": 85, "duration": 30, "unit": "min", "type": "hotfill"},
                {"temperature_C": 25, "duration": 6, "unit": "months", "type": "storage"},
            ],
        },
        {
            "code": "setoff_storage",
            "name": "Set-off + Storage",
            "description": "Manufacturing set-off followed by consumer use",
            "steps": [
                {"temperature_C": 25, "duration": 3, "unit": "months", "with_food": False},
                {"temperature_C": 4, "duration": 6, "unit": "months", "with_food": True},
            ],
        },
    ]
    return JSONResponse({
        "success": True,
        "presets": presets,
    })
