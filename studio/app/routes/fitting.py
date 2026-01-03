"""
Fitting Routes - Parameter Fitting for D and k

API endpoints for fitting diffusivity (D) and partition coefficient (k)
from experimental or synthetic migration data.

Based on example4.py workflow:
1. Define layer and food parameters
2. Generate synthetic data (with noise) OR load real experimental data
3. Fit D and k using optimization
4. Evaluate fit quality and compare with true values

@author: SFPPy Studio
@license: MIT
"""

import sys
from pathlib import Path
from typing import Optional, List, Dict, Any
from datetime import datetime
import json
import uuid

import numpy as np
from fastapi import APIRouter, HTTPException, Query, UploadFile, File
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

# Add parent paths
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

router = APIRouter()


# =============================================================================
# Data Models
# =============================================================================

class FittingLayerConfig(BaseModel):
    """Layer configuration for fitting"""
    polymer: str = Field(default="LDPE", description="Polymer type")
    thickness_um: float = Field(default=100, gt=0, description="Thickness in µm")
    C0: float = Field(default=1000, ge=0, description="Initial concentration (mg/kg)")
    D_true: Optional[float] = Field(None, description="True D for synthetic data (m²/s)")
    k_true: Optional[float] = Field(None, description="True k for synthetic data")


class FittingFoodConfig(BaseModel):
    """Food/medium configuration for fitting"""
    contact_time_days: float = Field(default=10, gt=0, description="Contact time in days")
    volume_L: float = Field(default=1.0, gt=0, description="Volume in liters")
    surface_area_dm2: Optional[float] = Field(None, gt=0, description="Surface area in dm²")
    surface_dm2: Optional[float] = Field(None, gt=0, description="Surface area in dm² (alias)")
    h_m_s: float = Field(default=1e-6, gt=0, description="Mass transfer coefficient (m/s)")
    CF0: float = Field(default=0, ge=0, description="Initial concentration in food (mg/kg)")
    k_food: float = Field(default=1.0, gt=0, description="Partition coefficient of food")

    @property
    def effective_surface_dm2(self) -> float:
        """Get surface area from either field name"""
        return self.surface_area_dm2 or self.surface_dm2 or 6.0


class SyntheticDataConfig(BaseModel):
    """Configuration for synthetic data generation"""
    layer: FittingLayerConfig = Field(default_factory=FittingLayerConfig)
    food: FittingFoodConfig = Field(default_factory=FittingFoodConfig)
    n_points: int = Field(default=30, ge=5, le=200, description="Number of data points")
    std_relative: Optional[float] = Field(None, ge=0, le=0.5, description="Relative noise std")
    noise_std: Optional[float] = Field(None, ge=0, le=0.5, description="Relative noise std (alias)")
    contact_time_days: Optional[float] = Field(None, gt=0, description="Contact time in days (top-level)")
    seed: Optional[int] = Field(None, description="Random seed for reproducibility")

    @property
    def effective_std_relative(self) -> float:
        """Get noise std from either field name"""
        return self.std_relative or self.noise_std or 0.01

    @property
    def effective_contact_time_days(self) -> float:
        """Get contact time from either location"""
        return self.contact_time_days or self.food.contact_time_days or 10.0


class ExperimentalDataPoint(BaseModel):
    """Single experimental data point"""
    time_days: float = Field(..., ge=0)
    CF: float = Field(..., ge=0, description="Concentration in food (mg/kg)")
    CF_std: Optional[float] = Field(None, ge=0, description="Standard deviation")


class ExperimentalDataInput(BaseModel):
    """Experimental data for fitting"""
    layer: FittingLayerConfig = Field(default_factory=FittingLayerConfig)
    food: FittingFoodConfig = Field(default_factory=FittingFoodConfig)
    data_points: List[ExperimentalDataPoint] = Field(..., min_length=3)


class FittingBounds(BaseModel):
    """Parameter bounds for fitting"""
    D_min: float = Field(default=1e-18, gt=0, description="Minimum D (m²/s)")
    D_max: float = Field(default=1e-10, gt=0, description="Maximum D (m²/s)")
    k_min: float = Field(default=0.001, gt=0, description="Minimum k")
    k_max: float = Field(default=1000, gt=0, description="Maximum k")


class FitRequest(BaseModel):
    """Request to perform fitting"""
    job_id: str = Field(..., description="Fitting job ID (from synthetic or real data)")
    bounds: FittingBounds = Field(default_factory=FittingBounds)
    fit_D: bool = Field(default=True, description="Fit diffusivity D")
    fit_k: bool = Field(default=True, description="Fit partition coefficient k")
    method: str = Field(default="L-BFGS-B", description="Optimization method")
    max_iterations: int = Field(default=1000, ge=10, le=10000)


# =============================================================================
# In-memory storage for fitting jobs
# =============================================================================

_fitting_jobs: Dict[str, Dict[str, Any]] = {}


# =============================================================================
# Helper Functions
# =============================================================================

def _generate_synthetic_migration_data(
    D_true: float,  # m²/s
    k_true: float,
    layer_thickness_m: float,
    C0: float,
    contact_time_s: float,
    volume_m3: float,
    surface_area_m2: float,
    n_points: int = 30,
    std_relative: float = 0.01,
    seed: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Generate synthetic migration kinetics data.

    Uses simplified analytical solution for short times or numerical approximation.
    Returns time points and CF values with added Gaussian noise.
    """
    if seed is not None:
        np.random.seed(seed)

    # Time points (logarithmically spaced for better coverage)
    t_min = contact_time_s / 1000
    t_max = contact_time_s
    t_s = np.logspace(np.log10(t_min), np.log10(t_max), n_points)

    # Equilibrium concentration
    # CF_eq = C0 * VP / (VP + k * VF)
    VP = layer_thickness_m * surface_area_m2  # Polymer volume approximation
    VF = volume_m3  # Food volume
    CF_eq = C0 * VP / (VP + k_true * VF)

    # Simplified kinetic model (exponential approach to equilibrium)
    # tau = l^2 / (D * pi^2) for diffusion-limited case
    tau = layer_thickness_m**2 / (D_true * np.pi**2)

    # CF(t) = CF_eq * (1 - exp(-t/tau))
    CF_clean = CF_eq * (1 - np.exp(-t_s / tau))

    # Add Gaussian noise
    noise = np.random.normal(0, std_relative * CF_eq, n_points)
    CF_noisy = np.maximum(0, CF_clean + noise)

    return {
        "time_s": t_s.tolist(),
        "time_days": (t_s / 86400).tolist(),
        "CF_clean": CF_clean.tolist(),
        "CF_noisy": CF_noisy.tolist(),
        "CF_equilibrium": float(CF_eq),
        "tau_s": float(tau),
        "D_true": D_true,
        "k_true": k_true,
    }


def _compute_model_CF(
    t_s: np.ndarray,
    D: float,
    k: float,
    layer_thickness_m: float,
    C0: float,
    volume_m3: float,
    surface_area_m2: float,
) -> np.ndarray:
    """Compute model CF values for given D and k"""
    VP = layer_thickness_m * surface_area_m2
    VF = volume_m3
    CF_eq = C0 * VP / (VP + k * VF)
    tau = layer_thickness_m**2 / (D * np.pi**2)
    CF = CF_eq * (1 - np.exp(-t_s / tau))
    return CF


def _objective_function(
    params: np.ndarray,
    t_s: np.ndarray,
    CF_exp: np.ndarray,
    layer_thickness_m: float,
    C0: float,
    volume_m3: float,
    surface_area_m2: float,
    fit_D: bool,
    fit_k: bool,
    D_fixed: Optional[float] = None,
    k_fixed: Optional[float] = None,
) -> float:
    """Objective function: sum of squared residuals"""
    idx = 0
    if fit_D:
        D = 10**params[idx]  # D is fit in log scale
        idx += 1
    else:
        D = D_fixed

    if fit_k:
        k = 10**params[idx]  # k is fit in log scale
    else:
        k = k_fixed

    CF_model = _compute_model_CF(
        t_s, D, k, layer_thickness_m, C0, volume_m3, surface_area_m2
    )

    residuals = CF_exp - CF_model
    return np.sum(residuals**2)


def _perform_fitting(
    t_s: np.ndarray,
    CF_exp: np.ndarray,
    layer_thickness_m: float,
    C0: float,
    volume_m3: float,
    surface_area_m2: float,
    bounds: FittingBounds,
    fit_D: bool = True,
    fit_k: bool = True,
    D_initial: Optional[float] = None,
    k_initial: Optional[float] = None,
    method: str = "L-BFGS-B",
    max_iterations: int = 1000,
) -> Dict[str, Any]:
    """Perform parameter fitting using scipy.optimize"""
    from scipy.optimize import minimize

    # Initial guesses (log scale)
    x0 = []
    bounds_list = []

    if fit_D:
        D_init = D_initial or 1e-14
        x0.append(np.log10(D_init))
        bounds_list.append((np.log10(bounds.D_min), np.log10(bounds.D_max)))

    if fit_k:
        k_init = k_initial or 0.1
        x0.append(np.log10(k_init))
        bounds_list.append((np.log10(bounds.k_min), np.log10(bounds.k_max)))

    x0 = np.array(x0)

    # Perform optimization
    result = minimize(
        _objective_function,
        x0,
        args=(
            t_s, CF_exp, layer_thickness_m, C0, volume_m3, surface_area_m2,
            fit_D, fit_k, D_initial if not fit_D else None, k_initial if not fit_k else None
        ),
        method=method,
        bounds=bounds_list,
        options={"maxiter": max_iterations},
    )

    # Extract fitted parameters
    idx = 0
    D_fitted = None
    k_fitted = None

    if fit_D:
        D_fitted = 10**result.x[idx]
        idx += 1

    if fit_k:
        k_fitted = 10**result.x[idx]

    # Compute model with fitted parameters
    D_final = D_fitted if fit_D else D_initial
    k_final = k_fitted if fit_k else k_initial

    CF_fitted = _compute_model_CF(
        t_s, D_final, k_final, layer_thickness_m, C0, volume_m3, surface_area_m2
    )

    # Statistics
    residuals = CF_exp - CF_fitted
    SS_res = np.sum(residuals**2)
    SS_tot = np.sum((CF_exp - np.mean(CF_exp))**2)
    R2 = 1 - SS_res / SS_tot if SS_tot > 0 else 0

    return {
        "D_fitted": D_fitted,
        "k_fitted": k_fitted,
        "D_fitted_cm2_s": D_fitted * 1e4 if D_fitted else None,  # Convert to cm²/s
        "k_fitted": k_fitted,
        "CF_fitted": CF_fitted.tolist(),
        "residuals": residuals.tolist(),
        "SS_res": float(SS_res),
        "R2": float(R2),
        "RMSE": float(np.sqrt(SS_res / len(CF_exp))),
        "optimization": {
            "success": result.success,
            "message": result.message,
            "n_iterations": result.nit if hasattr(result, 'nit') else None,
            "n_function_evals": result.nfev if hasattr(result, 'nfev') else None,
        },
    }


# =============================================================================
# Endpoints
# =============================================================================

@router.post("/synthetic/generate")
async def generate_synthetic_data(config: SyntheticDataConfig):
    """
    Generate synthetic migration data with prescribed D and k values.

    This simulates experimental data for testing fitting algorithms.
    Returns time points and CF values with added Gaussian noise.
    """
    # Validate D_true and k_true are provided
    if config.layer.D_true is None:
        config.layer.D_true = 1e-14  # Default: 1e-14 m²/s = 1e-10 cm²/s

    if config.layer.k_true is None:
        config.layer.k_true = 0.1  # Default k

    # Convert units - use effective properties for aliased fields
    layer_thickness_m = config.layer.thickness_um * 1e-6
    contact_time_days = config.effective_contact_time_days
    contact_time_s = contact_time_days * 86400
    volume_m3 = config.food.volume_L * 1e-3
    surface_area_m2 = config.food.effective_surface_dm2 * 0.01

    # Generate data
    data = _generate_synthetic_migration_data(
        D_true=config.layer.D_true,
        k_true=config.layer.k_true,
        layer_thickness_m=layer_thickness_m,
        C0=config.layer.C0,
        contact_time_s=contact_time_s,
        volume_m3=volume_m3,
        surface_area_m2=surface_area_m2,
        n_points=config.n_points,
        std_relative=config.effective_std_relative,
        seed=config.seed,
    )

    # Create job ID
    job_id = str(uuid.uuid4())[:8]

    # Store job data
    _fitting_jobs[job_id] = {
        "type": "synthetic",
        "created_at": datetime.utcnow().isoformat(),
        "config": config.model_dump(),
        "data": data,
        "layer_thickness_m": layer_thickness_m,
        "volume_m3": volume_m3,
        "surface_area_m2": surface_area_m2,
        "fit_result": None,
    }

    return JSONResponse({
        "success": True,
        "job_id": job_id,
        "config": {
            "D_true_m2_s": config.layer.D_true,
            "D_true_cm2_s": config.layer.D_true * 1e4,
            "k_true": config.layer.k_true,
            "n_points": config.n_points,
            "std_relative": config.effective_std_relative,
        },
        # Flat structure for JavaScript compatibility
        "time_days": data["time_days"],
        "CF_clean": data["CF_clean"],
        "CF_noisy": data["CF_noisy"],
        "CF_equilibrium": data["CF_equilibrium"],
        # Also nested for API consistency
        "data": {
            "time_days": data["time_days"],
            "CF_clean": data["CF_clean"],
            "CF_noisy": data["CF_noisy"],
            "CF_equilibrium": data["CF_equilibrium"],
        },
        "message": f"Synthetic data generated with {config.n_points} points",
    })


@router.post("/experimental/load")
async def load_experimental_data(data: ExperimentalDataInput):
    """
    Load experimental migration data for fitting.

    Accepts time-CF pairs from real experiments.
    """
    # Convert units
    layer_thickness_m = data.layer.thickness_um * 1e-6
    contact_time_s = data.food.contact_time_days * 86400
    volume_m3 = data.food.volume_L * 1e-3
    surface_area_m2 = data.food.effective_surface_dm2 * 0.01

    # Extract data points
    time_days = [p.time_days for p in data.data_points]
    time_s = [t * 86400 for t in time_days]
    CF_values = [p.CF for p in data.data_points]
    CF_std = [p.CF_std for p in data.data_points if p.CF_std is not None]

    # Create job ID
    job_id = str(uuid.uuid4())[:8]

    # Store job data
    _fitting_jobs[job_id] = {
        "type": "experimental",
        "created_at": datetime.utcnow().isoformat(),
        "config": {
            "layer": data.layer.model_dump(),
            "food": data.food.model_dump(),
        },
        "data": {
            "time_s": time_s,
            "time_days": time_days,
            "CF_exp": CF_values,
            "CF_std": CF_std if CF_std else None,
        },
        "layer_thickness_m": layer_thickness_m,
        "volume_m3": volume_m3,
        "surface_area_m2": surface_area_m2,
        "fit_result": None,
    }

    return JSONResponse({
        "success": True,
        "job_id": job_id,
        "n_points": len(data.data_points),
        "time_range_days": [min(time_days), max(time_days)],
        "CF_range": [min(CF_values), max(CF_values)],
        # Flat structure for JavaScript
        "time_days": time_days,
        "CF": CF_values,
        "message": f"Loaded {len(data.data_points)} experimental data points",
    })


class CSVDataInput(BaseModel):
    """Simple CSV data input for experimental data"""
    data_csv: str = Field(..., description="CSV data as string (time_days,CF per line)")
    thickness_um: float = Field(default=100, gt=0, description="Layer thickness in µm")
    C0: float = Field(default=1000, ge=0, description="Initial concentration (mg/kg)")
    volume_L: float = Field(default=1.0, gt=0, description="Volume in liters")
    surface_dm2: float = Field(default=6.0, gt=0, description="Surface area in dm²")


@router.post("/experimental/csv")
async def load_csv_data(data: CSVDataInput):
    """
    Load experimental data from CSV string.

    Expected format:
    time_days,CF
    0.5,12.3
    1.0,24.5
    ...
    """
    import csv
    from io import StringIO

    data_points = []
    reader = csv.reader(StringIO(data.data_csv.strip()))

    for row in reader:
        if len(row) >= 2:
            try:
                time_val = float(row[0].strip())
                cf_val = float(row[1].strip())
                # Skip header row
                if time_val >= 0 and cf_val >= 0:
                    data_points.append({"time_days": time_val, "CF": cf_val})
            except ValueError:
                continue

    if len(data_points) < 3:
        raise HTTPException(status_code=400, detail="Need at least 3 valid data points")

    time_days = [p["time_days"] for p in data_points]
    time_s = [t * 86400 for t in time_days]
    CF_values = [p["CF"] for p in data_points]

    # Create job
    job_id = str(uuid.uuid4())[:8]
    layer_thickness_m = data.thickness_um * 1e-6
    volume_m3 = data.volume_L * 1e-3
    surface_area_m2 = data.surface_dm2 * 0.01

    _fitting_jobs[job_id] = {
        "type": "experimental",
        "created_at": datetime.utcnow().isoformat(),
        "config": {
            "layer": {"thickness_um": data.thickness_um, "C0": data.C0},
            "food": {"volume_L": data.volume_L, "surface_dm2": data.surface_dm2},
        },
        "data": {
            "time_s": time_s,
            "time_days": time_days,
            "CF_exp": CF_values,
        },
        "layer_thickness_m": layer_thickness_m,
        "volume_m3": volume_m3,
        "surface_area_m2": surface_area_m2,
        "fit_result": None,
    }

    return JSONResponse({
        "success": True,
        "job_id": job_id,
        "n_points": len(data_points),
        "time_days": time_days,
        "CF": CF_values,
        "message": f"Loaded {len(data_points)} data points from CSV",
    })


@router.post("/experimental/upload")
async def upload_experimental_data(
    file: UploadFile = File(...),
    thickness_um: float = Query(100, gt=0),
    C0: float = Query(1000, ge=0),
    volume_L: float = Query(1.0, gt=0),
    surface_area_dm2: float = Query(6.0, gt=0),
):
    """
    Upload experimental data from CSV file.

    Expected CSV format:
    time_days,CF
    0.5,12.3
    1.0,24.5
    ...
    """
    import csv
    from io import StringIO

    content = await file.read()
    text = content.decode("utf-8")

    data_points = []
    reader = csv.DictReader(StringIO(text))

    for row in reader:
        try:
            time_days = float(row.get("time_days") or row.get("time"))
            CF = float(row.get("CF") or row.get("concentration"))
            CF_std = float(row["CF_std"]) if "CF_std" in row and row["CF_std"] else None
            data_points.append({"time_days": time_days, "CF": CF, "CF_std": CF_std})
        except (ValueError, KeyError) as e:
            continue

    if len(data_points) < 3:
        raise HTTPException(status_code=400, detail="Need at least 3 valid data points")

    # Convert to ExperimentalDataInput and process
    layer_config = FittingLayerConfig(thickness_um=thickness_um, C0=C0)
    food_config = FittingFoodConfig(
        volume_L=volume_L,
        surface_area_dm2=surface_area_dm2,
        contact_time_days=max(p["time_days"] for p in data_points),
    )

    # Create job
    job_id = str(uuid.uuid4())[:8]
    layer_thickness_m = thickness_um * 1e-6
    volume_m3 = volume_L * 1e-3
    surface_area_m2 = surface_area_dm2 * 0.01

    _fitting_jobs[job_id] = {
        "type": "experimental",
        "created_at": datetime.utcnow().isoformat(),
        "config": {
            "layer": layer_config.model_dump(),
            "food": food_config.model_dump(),
            "source_file": file.filename,
        },
        "data": {
            "time_s": [p["time_days"] * 86400 for p in data_points],
            "time_days": [p["time_days"] for p in data_points],
            "CF_exp": [p["CF"] for p in data_points],
            "CF_std": [p.get("CF_std") for p in data_points],
        },
        "layer_thickness_m": layer_thickness_m,
        "volume_m3": volume_m3,
        "surface_area_m2": surface_area_m2,
        "fit_result": None,
    }

    return JSONResponse({
        "success": True,
        "job_id": job_id,
        "filename": file.filename,
        "n_points": len(data_points),
        "message": f"Loaded {len(data_points)} points from {file.filename}",
    })


@router.post("/fit")
async def perform_fit(request: FitRequest):
    """
    Perform D and k fitting on loaded data.

    Uses scipy.optimize to minimize the squared error between
    model predictions and experimental/synthetic data.
    """
    job_id = request.job_id

    if job_id not in _fitting_jobs:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    job = _fitting_jobs[job_id]

    # Get data
    if job["type"] == "synthetic":
        t_s = np.array(job["data"]["time_s"])
        CF_exp = np.array(job["data"]["CF_noisy"])
        C0 = job["config"]["layer"]["C0"]
    else:
        t_s = np.array(job["data"]["time_s"])
        CF_exp = np.array(job["data"]["CF_exp"])
        C0 = job["config"]["layer"]["C0"]

    # Perform fitting
    try:
        result = _perform_fitting(
            t_s=t_s,
            CF_exp=CF_exp,
            layer_thickness_m=job["layer_thickness_m"],
            C0=C0,
            volume_m3=job["volume_m3"],
            surface_area_m2=job["surface_area_m2"],
            bounds=request.bounds,
            fit_D=request.fit_D,
            fit_k=request.fit_k,
            method=request.method,
            max_iterations=request.max_iterations,
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Fitting failed: {str(e)}")

    # Store result
    job["fit_result"] = result

    # Compare with true values if synthetic
    comparison = None
    if job["type"] == "synthetic" and job["data"].get("D_true"):
        D_true = job["data"]["D_true"]
        k_true = job["data"]["k_true"]
        comparison = {
            "D_error_percent": 100 * abs(result["D_fitted"] - D_true) / D_true if result["D_fitted"] else None,
            "k_error_percent": 100 * abs(result["k_fitted"] - k_true) / k_true if result["k_fitted"] else None,
            "D_true": D_true,
            "k_true": k_true,
        }

    # Build fitted curve with time data
    time_days = job["data"]["time_days"]

    return JSONResponse({
        "success": True,
        "job_id": job_id,
        # Flat structure expected by JavaScript
        "D_fitted": result["D_fitted"],
        "k_fitted": result["k_fitted"],
        "R2": result["R2"],
        "RMSE": result["RMSE"],
        # Nested structures for completeness
        "fitted_parameters": {
            "D_m2_s": result["D_fitted"],
            "D_cm2_s": result["D_fitted_cm2_s"],
            "k": result["k_fitted"],
        },
        "fit_quality": {
            "R2": result["R2"],
            "RMSE": result["RMSE"],
            "SS_res": result["SS_res"],
        },
        "comparison": comparison,  # JavaScript expects 'comparison' not 'comparison_with_true'
        "comparison_with_true": comparison,
        "optimization": result["optimization"],
        "fitted_curve": {
            "time_days": time_days,
            "CF_fitted": result["CF_fitted"],
        },
        "CF_fitted": result["CF_fitted"],
    })


@router.get("/job/{job_id}")
async def get_fitting_job(job_id: str):
    """Get details of a fitting job"""
    if job_id not in _fitting_jobs:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    job = _fitting_jobs[job_id]

    return JSONResponse({
        "success": True,
        "job_id": job_id,
        "type": job["type"],
        "created_at": job["created_at"],
        "config": job["config"],
        "data": {
            "time_days": job["data"]["time_days"],
            "CF": job["data"].get("CF_noisy") or job["data"].get("CF_exp"),
        },
        "has_fit_result": job["fit_result"] is not None,
        "fit_result": job["fit_result"],
    })


@router.get("/jobs")
async def list_fitting_jobs():
    """List all fitting jobs"""
    jobs = []
    for job_id, job in _fitting_jobs.items():
        jobs.append({
            "job_id": job_id,
            "type": job["type"],
            "created_at": job["created_at"],
            "has_fit_result": job["fit_result"] is not None,
        })

    return JSONResponse({
        "success": True,
        "jobs": jobs,
        "count": len(jobs),
    })


@router.delete("/job/{job_id}")
async def delete_fitting_job(job_id: str):
    """Delete a fitting job"""
    if job_id not in _fitting_jobs:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    del _fitting_jobs[job_id]

    return JSONResponse({
        "success": True,
        "message": f"Job {job_id} deleted",
    })


@router.get("/presets")
async def get_fitting_presets():
    """Get predefined fitting presets (from example4.py)"""
    return JSONResponse({
        "success": True,
        "presets": [
            {
                "name": "Example 4 Default",
                "description": "100 µm layer, D=1e-14 m²/s, k=0.1, 10 days contact",
                "code": "example4",
                "config": {
                    "layer": {
                        "polymer": "generic",
                        "thickness_um": 100,
                        "C0": 1000,
                        "D_true": 1e-14,
                        "k_true": 0.1,
                    },
                    "food": {
                        "contact_time_days": 10,
                        "volume_L": 1.0,
                        "surface_area_dm2": 6.0,
                        "h_m_s": 1e-6,
                        "CF0": 0,
                        "k_food": 1.0,
                    },
                    "n_points": 30,
                    "std_relative": 0.01,
                },
            },
            {
                "name": "High Noise",
                "description": "Same as Example 4 but with 5% noise",
                "code": "high_noise",
                "config": {
                    "layer": {
                        "polymer": "generic",
                        "thickness_um": 100,
                        "C0": 1000,
                        "D_true": 1e-14,
                        "k_true": 0.1,
                    },
                    "food": {
                        "contact_time_days": 10,
                        "volume_L": 1.0,
                        "surface_area_dm2": 6.0,
                        "h_m_s": 1e-6,
                        "CF0": 0,
                        "k_food": 1.0,
                    },
                    "n_points": 30,
                    "std_relative": 0.05,
                },
            },
            {
                "name": "Fast Diffusion",
                "description": "Higher D for fast migrating substances",
                "code": "fast_diffusion",
                "config": {
                    "layer": {
                        "polymer": "LDPE",
                        "thickness_um": 100,
                        "C0": 1000,
                        "D_true": 1e-12,
                        "k_true": 0.5,
                    },
                    "food": {
                        "contact_time_days": 2,
                        "volume_L": 0.5,
                        "surface_area_dm2": 4.0,
                        "h_m_s": 1e-6,
                        "CF0": 0,
                        "k_food": 1.0,
                    },
                    "n_points": 20,
                    "std_relative": 0.02,
                },
            },
            {
                "name": "Slow Diffusion",
                "description": "Lower D for slow migrating substances (barrier materials)",
                "code": "slow_diffusion",
                "config": {
                    "layer": {
                        "polymer": "PET",
                        "thickness_um": 50,
                        "C0": 500,
                        "D_true": 1e-16,
                        "k_true": 0.01,
                    },
                    "food": {
                        "contact_time_days": 30,
                        "volume_L": 1.0,
                        "surface_area_dm2": 6.0,
                        "h_m_s": 1e-6,
                        "CF0": 0,
                        "k_food": 1.0,
                    },
                    "n_points": 40,
                    "std_relative": 0.01,
                },
            },
        ],
    })
