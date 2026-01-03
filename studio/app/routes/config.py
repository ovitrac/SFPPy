"""
Configuration Routes - Application Settings

API endpoints for user preferences, solver settings, and application configuration.
"""

import sys
import json
import socket
import platform
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, Any

from fastapi import APIRouter, HTTPException
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

# Add parent paths
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

router = APIRouter()

# Config directory
CONFIG_DIR = Path(__file__).parent.parent.parent / "config"
CONFIG_DIR.mkdir(exist_ok=True)

USER_CONFIG_FILE = CONFIG_DIR / "user_config.json"


# ========== DATA MODELS ==========

class UserProfile(BaseModel):
    name: str = Field(default="", description="User name for traceability")
    organization: str = Field(default="", description="Organization name")
    email: str = Field(default="", description="Contact email")
    role: str = Field(default="analyst", description="User role")


class SolverSettings(BaseModel):
    """
    Solver settings aligned with patankar.useroverride parameters.

    These settings control the Patankar finite-volume solver behavior:
    - ntimes: Number of stored simulation time points (affects resolution of CF(t))
    - nmesh: Total FV nodes across all layers (affects spatial resolution)
    - nmeshmin: Minimum FV nodes per layer (ensures thin layers are resolved)
    - timescale: Time point distribution ("sqrt" for diffusion, "linear" for uniform)
    - RelTol/AbsTol: ODE integration tolerances
    - nmax: Number of concentration profiles to store for C(x) visualization
    """
    # Patankar migration settings (from useroverride.py)
    ntimes: int = Field(default=1000, ge=50, le=20000, description="Number of stored simulation times")
    nmesh: int = Field(default=600, ge=100, le=5000, description="Total number of FV volumes in assembly")
    nmeshmin: int = Field(default=20, ge=5, le=100, description="Minimum FV volumes per layer")
    timescale: str = Field(default="sqrt", description="Time distribution: sqrt (diffusion) or linear")
    RelTol: float = Field(default=1e-6, ge=1e-12, le=1e-2, description="Relative tolerance for PDE integration")
    AbsTol: float = Field(default=1e-6, ge=1e-12, le=1e-2, description="Absolute tolerance for PDE integration")
    nmax: int = Field(default=15, ge=5, le=50, description="Number of concentration profiles per substance")

    # Legacy compatibility fields (mapped to new names)
    nx_default: int = Field(default=600, ge=10, le=5000, description="Alias for nmesh")
    nt_default: int = Field(default=1000, ge=100, le=20000, description="Alias for ntimes")
    atol: float = Field(default=1e-6, description="Alias for AbsTol")
    rtol: float = Field(default=1e-6, description="Alias for RelTol")
    max_iterations: int = Field(default=1000, description="Maximum solver iterations (not used in Patankar)")


class EstimatorSettings(BaseModel):
    D_model: str = Field(default="piringer", description="Diffusivity model: piringer, welle, measured")
    k_model: str = Field(default="fhp", description="Partition model: fhp, nernst, measured")
    temperature_correction: bool = Field(default=True, description="Apply Arrhenius correction")
    use_plasticization: bool = Field(default=True, description="Account for plasticizer effects")


class DisplaySettings(BaseModel):
    theme: str = Field(default="light", description="UI theme: light, dark, auto")
    units_thickness: str = Field(default="um", description="Default thickness unit")
    units_concentration: str = Field(default="mg/kg", description="Default concentration unit")
    units_time: str = Field(default="days", description="Default time unit")
    decimal_places: int = Field(default=3, description="Decimal places for display")
    scientific_notation_threshold: float = Field(default=1e-4, description="Switch to scientific notation below this")


class AppConfig(BaseModel):
    user: UserProfile = Field(default_factory=UserProfile)
    solver: SolverSettings = Field(default_factory=SolverSettings)
    estimators: EstimatorSettings = Field(default_factory=EstimatorSettings)
    display: DisplaySettings = Field(default_factory=DisplaySettings)


# ========== HELPER FUNCTIONS ==========

def load_config() -> AppConfig:
    """Load configuration from disk."""
    if USER_CONFIG_FILE.exists():
        try:
            with open(USER_CONFIG_FILE, 'r') as f:
                data = json.load(f)
            return AppConfig(**data)
        except (json.JSONDecodeError, TypeError):
            pass
    return AppConfig()


def save_config(config: AppConfig):
    """Save configuration to disk."""
    with open(USER_CONFIG_FILE, 'w') as f:
        json.dump(config.dict(), f, indent=2)


# ========== API ENDPOINTS ==========

@router.get("/")
async def get_config():
    """Get current configuration."""
    config = load_config()
    return JSONResponse({
        "success": True,
        "config": config.dict(),
    })


@router.put("/")
async def update_config(config: AppConfig):
    """Update full configuration."""
    save_config(config)
    return JSONResponse({
        "success": True,
        "message": "Configuration saved",
        "config": config.dict(),
    })


@router.get("/user")
async def get_user_profile():
    """Get user profile."""
    config = load_config()
    return JSONResponse({
        "success": True,
        "user": config.user.dict(),
    })


@router.put("/user")
async def update_user_profile(user: UserProfile):
    """Update user profile."""
    config = load_config()
    config.user = user
    save_config(config)
    return JSONResponse({
        "success": True,
        "message": "User profile updated",
        "user": user.dict(),
    })


@router.get("/solver")
async def get_solver_settings():
    """Get solver settings."""
    config = load_config()
    return JSONResponse({
        "success": True,
        "solver": config.solver.dict(),
    })


@router.put("/solver")
async def update_solver_settings(solver: SolverSettings):
    """Update solver settings."""
    config = load_config()
    config.solver = solver
    save_config(config)
    return JSONResponse({
        "success": True,
        "message": "Solver settings updated",
        "solver": solver.dict(),
    })


@router.get("/estimators")
async def get_estimator_settings():
    """Get estimator settings."""
    config = load_config()
    return JSONResponse({
        "success": True,
        "estimators": config.estimators.dict(),
    })


@router.put("/estimators")
async def update_estimator_settings(estimators: EstimatorSettings):
    """Update estimator settings."""
    config = load_config()
    config.estimators = estimators
    save_config(config)
    return JSONResponse({
        "success": True,
        "message": "Estimator settings updated",
        "estimators": estimators.dict(),
    })


@router.get("/display")
async def get_display_settings():
    """Get display settings."""
    config = load_config()
    return JSONResponse({
        "success": True,
        "display": config.display.dict(),
    })


@router.put("/display")
async def update_display_settings(display: DisplaySettings):
    """Update display settings."""
    config = load_config()
    config.display = display
    save_config(config)
    return JSONResponse({
        "success": True,
        "message": "Display settings updated",
        "display": display.dict(),
    })


@router.post("/reset")
async def reset_config(section: Optional[str] = None):
    """Reset configuration to defaults."""
    if section:
        config = load_config()
        if section == "user":
            config.user = UserProfile()
        elif section == "solver":
            config.solver = SolverSettings()
        elif section == "estimators":
            config.estimators = EstimatorSettings()
        elif section == "display":
            config.display = DisplaySettings()
        else:
            raise HTTPException(status_code=400, detail=f"Unknown section: {section}")
        save_config(config)
    else:
        config = AppConfig()
        save_config(config)

    return JSONResponse({
        "success": True,
        "message": f"Configuration {'section ' + section + ' ' if section else ''}reset to defaults",
        "config": config.dict(),
    })


@router.get("/traceability")
async def get_traceability_info():
    """Get full traceability information for reports."""
    config = load_config()

    # Get SFPPy version
    try:
        from patankar.migration import __version__ as sfppy_version
    except ImportError:
        sfppy_version = "unknown"

    # Get Studio version
    try:
        from studio.version import __version__ as studio_version
    except ImportError:
        studio_version = "unknown"

    return JSONResponse({
        "success": True,
        "traceability": {
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "user": config.user.dict(),
            "machine": {
                "hostname": socket.gethostname(),
                "platform": platform.system(),
                "platform_version": platform.version(),
                "architecture": platform.machine(),
                "python_version": platform.python_version(),
            },
            "software": {
                "sfppy_version": sfppy_version,
                "studio_version": studio_version,
            },
        },
    })


@router.get("/models/diffusivity")
async def list_diffusivity_models():
    """List available diffusivity estimation models."""
    models = [
        {
            "code": "piringer",
            "name": "Piringer Model",
            "description": "EU standard model for polyolefins and PET",
            "parameters": ["Ap", "tau"],
            "references": ["Piringer et al., Food Additives & Contaminants, 2008"],
            "default": True,
        },
        {
            "code": "welle",
            "name": "Welle Model",
            "description": "Extended model for various polymers",
            "parameters": ["A", "B", "C"],
            "references": ["Welle, Packaging Technology and Science, 2013"],
            "default": False,
        },
        {
            "code": "measured",
            "name": "Measured Values",
            "description": "Use experimentally determined D values",
            "parameters": [],
            "references": [],
            "default": False,
        },
    ]

    return JSONResponse({
        "success": True,
        "models": models,
    })


@router.get("/models/partition")
async def list_partition_models():
    """List available partition coefficient models."""
    models = [
        {
            "code": "fhp",
            "name": "Flory-Huggins-Prigogine",
            "description": "Thermodynamic model based on solubility parameters",
            "parameters": ["chi", "delta"],
            "references": ["Flory, J. Chem. Phys., 1942"],
            "default": True,
        },
        {
            "code": "nernst",
            "name": "Nernst Distribution",
            "description": "Simple equilibrium distribution",
            "parameters": ["K"],
            "references": [],
            "default": False,
        },
        {
            "code": "measured",
            "name": "Measured Values",
            "description": "Use experimentally determined k values",
            "parameters": [],
            "references": [],
            "default": False,
        },
    ]

    return JSONResponse({
        "success": True,
        "models": models,
    })


@router.get("/units")
async def list_available_units():
    """List all available units for conversion."""
    units = {
        "thickness": [
            {"code": "nm", "name": "Nanometers", "factor": 1e-9},
            {"code": "um", "name": "Micrometers", "factor": 1e-6},
            {"code": "mm", "name": "Millimeters", "factor": 1e-3},
            {"code": "cm", "name": "Centimeters", "factor": 1e-2},
            {"code": "m", "name": "Meters", "factor": 1.0},
        ],
        "concentration": [
            {"code": "mg/kg", "name": "mg per kg", "factor": 1e-6},
            {"code": "ppm", "name": "Parts per million", "factor": 1e-6},
            {"code": "g/kg", "name": "g per kg", "factor": 1e-3},
            {"code": "%", "name": "Percent", "factor": 1e-2},
        ],
        "time": [
            {"code": "s", "name": "Seconds", "factor": 1},
            {"code": "min", "name": "Minutes", "factor": 60},
            {"code": "h", "name": "Hours", "factor": 3600},
            {"code": "days", "name": "Days", "factor": 86400},
            {"code": "weeks", "name": "Weeks", "factor": 604800},
            {"code": "months", "name": "Months (30d)", "factor": 2592000},
            {"code": "years", "name": "Years (365d)", "factor": 31536000},
        ],
        "diffusivity": [
            {"code": "m2/s", "name": "m²/s", "factor": 1.0},
            {"code": "cm2/s", "name": "cm²/s", "factor": 1e-4},
        ],
        "temperature": [
            {"code": "C", "name": "Celsius", "offset": 273.15},
            {"code": "K", "name": "Kelvin", "offset": 0},
            {"code": "F", "name": "Fahrenheit", "formula": "(T-32)*5/9+273.15"},
        ],
    }

    return JSONResponse({
        "success": True,
        "units": units,
    })


@router.get("/regulations")
async def list_regulations():
    """List supported regulatory frameworks."""
    regulations = [
        {
            "code": "EU_10_2011",
            "name": "EU 10/2011",
            "description": "European regulation on plastic materials for food contact",
            "region": "European Union",
            "default_SML": 60.0,
            "default_OML": 10.0,
        },
        {
            "code": "FDA_21CFR",
            "name": "FDA 21 CFR",
            "description": "US FDA regulations for food contact substances",
            "region": "United States",
            "default_SML": None,
        },
        {
            "code": "GB_4806",
            "name": "GB 4806 Series",
            "description": "Chinese national standards for food contact materials",
            "region": "China",
            "default_SML": None,
        },
    ]

    return JSONResponse({
        "success": True,
        "regulations": regulations,
    })
