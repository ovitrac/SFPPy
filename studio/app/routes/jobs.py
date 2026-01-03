"""
Jobs Routes - Simulation History Management

API endpoints for listing, loading, and managing past simulations.
"""

import sys
import json
import shutil
from pathlib import Path
from datetime import datetime
from typing import List, Optional

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

# Add parent paths
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

router = APIRouter()

# Jobs directory
JOBS_DIR = Path(__file__).parent.parent.parent / "jobs"
JOBS_DIR.mkdir(exist_ok=True)

EXPORTS_DIR = Path(__file__).parent.parent.parent / "exports"
EXPORTS_DIR.mkdir(exist_ok=True)


# ========== DATA MODELS ==========

class JobSummary(BaseModel):
    job_id: str
    name: str
    status: str
    created_at: str
    completed_at: Optional[str]
    layer_count: int
    step_count: int


class ExportFormat(BaseModel):
    format: str = Field(..., description="Export format: json, csv, xlsx, pdf")
    include_config: bool = True
    include_results: bool = True
    include_plots: bool = False


# ========== HELPER FUNCTIONS ==========

def list_job_files() -> List[Path]:
    """List all job files."""
    return sorted(JOBS_DIR.glob("*.json"), key=lambda p: p.stat().st_mtime, reverse=True)


def load_job_file(job_id: str) -> Optional[dict]:
    """Load a job file by ID."""
    job_file = JOBS_DIR / f"{job_id}.json"
    if job_file.exists():
        with open(job_file, 'r') as f:
            return json.load(f)
    return None


def delete_job_file(job_id: str) -> bool:
    """Delete a job file."""
    job_file = JOBS_DIR / f"{job_id}.json"
    if job_file.exists():
        job_file.unlink()
        return True
    return False


# ========== API ENDPOINTS ==========

@router.get("/list")
async def list_jobs(
    status: Optional[str] = Query(None, description="Filter by status"),
    limit: int = Query(50, ge=1, le=200),
    offset: int = Query(0, ge=0),
):
    """List all simulation jobs."""
    jobs = []

    for job_file in list_job_files():
        try:
            with open(job_file, 'r') as f:
                data = json.load(f)

            # Apply status filter
            if status and data.get("status") != status:
                continue

            # Extract summary
            config = data.get("input_config", {})
            jobs.append({
                "job_id": data.get("job_id"),
                "name": data.get("name", "Untitled"),
                "status": data.get("status"),
                "created_at": data.get("created_at"),
                "completed_at": data.get("completed_at"),
                "layer_count": len(config.get("layers", [])),
                "step_count": len(config.get("steps", [])),
                "has_results": data.get("results") is not None,
            })

        except (json.JSONDecodeError, KeyError):
            continue

    # Apply pagination
    total = len(jobs)
    jobs = jobs[offset:offset + limit]

    return JSONResponse({
        "success": True,
        "jobs": jobs,
        "total": total,
        "offset": offset,
        "limit": limit,
    })


@router.get("/detail/{job_id}")
async def get_job_detail(job_id: str):
    """Get full details of a job including config and results."""
    job = load_job_file(job_id)

    if not job:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    return JSONResponse({
        "success": True,
        **job,
    })


@router.delete("/delete/{job_id}")
async def delete_job(job_id: str):
    """Delete a job and its results."""
    if delete_job_file(job_id):
        return JSONResponse({
            "success": True,
            "message": f"Job {job_id} deleted",
        })

    raise HTTPException(status_code=404, detail=f"Job {job_id} not found")


@router.delete("/clear")
async def clear_old_jobs(
    older_than_days: int = Query(30, ge=1, description="Delete jobs older than N days"),
    status: Optional[str] = Query(None, description="Only delete jobs with this status"),
):
    """Clear old jobs from storage."""
    deleted = 0
    cutoff = datetime.utcnow().timestamp() - (older_than_days * 86400)

    for job_file in list_job_files():
        try:
            if job_file.stat().st_mtime < cutoff:
                if status:
                    with open(job_file, 'r') as f:
                        data = json.load(f)
                    if data.get("status") != status:
                        continue

                job_file.unlink()
                deleted += 1

        except (json.JSONDecodeError, PermissionError):
            continue

    return JSONResponse({
        "success": True,
        "deleted_count": deleted,
        "older_than_days": older_than_days,
    })


@router.post("/duplicate/{job_id}")
async def duplicate_job(job_id: str, new_name: Optional[str] = None):
    """Duplicate a job configuration for a new run."""
    import uuid

    job = load_job_file(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    # Create new job with same config
    new_job_id = str(uuid.uuid4())[:8]
    new_job = {
        "job_id": new_job_id,
        "name": new_name or f"Copy of {job.get('name', 'Untitled')}",
        "status": "draft",
        "created_at": datetime.utcnow().isoformat() + "Z",
        "started_at": None,
        "completed_at": None,
        "input_config": job.get("input_config"),
        "results": None,
        "error": None,
        "copied_from": job_id,
    }

    # Save new job
    new_job_file = JOBS_DIR / f"{new_job_id}.json"
    with open(new_job_file, 'w') as f:
        json.dump(new_job, f, indent=2)

    return JSONResponse({
        "success": True,
        "job_id": new_job_id,
        "copied_from": job_id,
        "message": "Job duplicated successfully",
    })


@router.post("/export/{job_id}")
async def export_job(job_id: str, export_format: str = "json"):
    """Export job results to various formats."""
    job = load_job_file(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    if job.get("status") != "completed":
        raise HTTPException(status_code=400, detail="Can only export completed jobs")

    timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    export_name = f"{job_id}_{timestamp}"

    if export_format == "json":
        export_file = EXPORTS_DIR / f"{export_name}.json"
        with open(export_file, 'w') as f:
            json.dump(job, f, indent=2)

        return JSONResponse({
            "success": True,
            "format": "json",
            "file": str(export_file.name),
            "path": str(export_file),
        })

    elif export_format == "csv":
        # Export results as CSV
        import csv

        results = job.get("results", {})
        export_file = EXPORTS_DIR / f"{export_name}.csv"

        with open(export_file, 'w', newline='') as f:
            writer = csv.writer(f)

            # Write time series data
            time_s = results.get("time_s", [])
            time_days = results.get("time_days", [])
            CF = results.get("CF_mg_kg", [])

            writer.writerow(["time_s", "time_days", "CF_mg_kg"])
            for i in range(len(time_s)):
                writer.writerow([
                    time_s[i] if i < len(time_s) else "",
                    time_days[i] if i < len(time_days) else "",
                    CF[i] if i < len(CF) else "",
                ])

        return JSONResponse({
            "success": True,
            "format": "csv",
            "file": str(export_file.name),
            "path": str(export_file),
        })

    else:
        raise HTTPException(status_code=400, detail=f"Unsupported export format: {export_format}")


@router.get("/compare")
async def compare_jobs(job_ids: str = Query(..., description="Comma-separated job IDs")):
    """Compare results from multiple jobs."""
    ids = [id.strip() for id in job_ids.split(",")]

    if len(ids) < 2:
        raise HTTPException(status_code=400, detail="Need at least 2 jobs to compare")

    if len(ids) > 5:
        raise HTTPException(status_code=400, detail="Maximum 5 jobs for comparison")

    comparison = []
    for job_id in ids:
        job = load_job_file(job_id)
        if not job:
            continue

        results = job.get("results", {})
        config = job.get("input_config", {})

        comparison.append({
            "job_id": job_id,
            "name": job.get("name"),
            "status": job.get("status"),
            "layer_count": len(config.get("layers", [])),
            "step_count": len(config.get("steps", [])),
            "final_CF": results.get("final_CF_mg_kg"),
            "compliant": results.get("compliant"),
            "SML": results.get("SML_mg_kg"),
            "is_mock": results.get("is_mock", False),
        })

    return JSONResponse({
        "success": True,
        "comparison": comparison,
        "count": len(comparison),
    })


@router.get("/stats")
async def get_job_stats():
    """Get statistics about stored jobs."""
    stats = {
        "total": 0,
        "by_status": {},
        "total_size_mb": 0,
        "oldest_job": None,
        "newest_job": None,
    }

    job_files = list_job_files()

    for job_file in job_files:
        try:
            stats["total_size_mb"] += job_file.stat().st_size / (1024 * 1024)

            with open(job_file, 'r') as f:
                data = json.load(f)

            status = data.get("status", "unknown")
            stats["by_status"][status] = stats["by_status"].get(status, 0) + 1
            stats["total"] += 1

            created = data.get("created_at")
            if created:
                if stats["oldest_job"] is None or created < stats["oldest_job"]:
                    stats["oldest_job"] = created
                if stats["newest_job"] is None or created > stats["newest_job"]:
                    stats["newest_job"] = created

        except (json.JSONDecodeError, KeyError):
            continue

    stats["total_size_mb"] = round(stats["total_size_mb"], 2)

    return JSONResponse({
        "success": True,
        **stats,
    })


@router.post("/save-draft")
async def save_draft(config: dict, name: str = "Draft"):
    """Save a configuration as a draft for later use."""
    import uuid

    job_id = str(uuid.uuid4())[:8]

    draft = {
        "job_id": job_id,
        "name": name,
        "status": "draft",
        "created_at": datetime.utcnow().isoformat() + "Z",
        "started_at": None,
        "completed_at": None,
        "input_config": config,
        "results": None,
        "error": None,
    }

    job_file = JOBS_DIR / f"{job_id}.json"
    with open(job_file, 'w') as f:
        json.dump(draft, f, indent=2)

    return JSONResponse({
        "success": True,
        "job_id": job_id,
        "message": "Draft saved successfully",
    })


@router.get("/load/{job_id}/config")
async def load_job_config(job_id: str):
    """Load only the configuration from a job (for re-running)."""
    job = load_job_file(job_id)

    if not job:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    return JSONResponse({
        "success": True,
        "job_id": job_id,
        "name": job.get("name"),
        "config": job.get("input_config"),
    })
