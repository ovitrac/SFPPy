"""
survey/app/main.py — FastAPI Application for Family Editor
==========================================================

Interactive web application for creating and editing substance families.
Features:
- PubChem substance lookup with molecule visualization
- Concentration distribution preview
- Family YAML import/export
- Validation and error reporting

Usage:
    uvicorn survey.app.main:app --reload --port 8000

@project: SFPPy/INSERM — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

import sys
import base64
import io
import traceback
from pathlib import Path
from typing import List, Dict, Any, Optional
from dataclasses import dataclass, field, asdict

import yaml
import numpy as np
from fastapi import FastAPI, Request, HTTPException, Form, UploadFile, File
from fastapi.responses import HTMLResponse, JSONResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from patankar.loadpubchem import migrant


# =============================================================================
# Data Models
# =============================================================================

@dataclass
class SubstanceInfo:
    """Substance information from PubChem."""
    name: str
    cas: Optional[str] = None
    mass_g_mol: float = 0.0
    formula: str = ""
    smiles: str = ""
    cid: Optional[int] = None
    logP: Optional[float] = None
    image_base64: str = ""
    C0_min: float = 0.0
    C0_likely: float = 100.0
    C0_max: float = 500.0
    notes: str = ""
    error: Optional[str] = None
    validated: bool = False

    def to_yaml_dict(self) -> Dict[str, Any]:
        """Convert to YAML-compatible dict."""
        d = {
            "name": self.name,
            "cas": self.cas,
            "full_name": self.formula,
            "C0_min": self.C0_min,
            "C0_likely": self.C0_likely,
            "C0_max": self.C0_max,
        }
        if self.notes:
            d["notes"] = self.notes
        return d


@dataclass
class FamilyInfo:
    """Family of substances."""
    name: str = "new_family"
    description: str = ""
    substances: List[SubstanceInfo] = field(default_factory=list)

    def to_yaml(self) -> str:
        """Convert to YAML string."""
        data = {
            "name": self.name,
            "description": self.description,
            "substances": [s.to_yaml_dict() for s in self.substances if s.validated],
        }
        return yaml.dump(data, default_flow_style=False, allow_unicode=True, sort_keys=False)


# =============================================================================
# In-memory storage (session-based)
# =============================================================================

# Simple in-memory storage for current family being edited
current_family: FamilyInfo = FamilyInfo()


# =============================================================================
# PubChem Lookup Functions
# =============================================================================

def lookup_substance(identifier: str) -> SubstanceInfo:
    """
    Look up substance by name or CAS number.

    Parameters
    ----------
    identifier : str
        Substance name or CAS number.

    Returns
    -------
    SubstanceInfo
        Populated substance info, with error field if lookup failed.
    """
    info = SubstanceInfo(name=identifier)

    try:
        # Try to fetch from PubChem
        m = migrant(name=identifier)

        info.name = m.compound or identifier
        info.mass_g_mol = float(m.M) if m.M else 0.0
        info.formula = m.formula or ""
        info.smiles = m.smiles or ""
        info.cid = m.cid

        # Get CAS (first one if multiple)
        if m.CAS:
            info.cas = m.CAS[0] if isinstance(m.CAS, list) else str(m.CAS)

        # Get logP
        if hasattr(m, 'logP') and m.logP is not None:
            try:
                info.logP = float(m.logP[0]) if hasattr(m.logP, '__iter__') else float(m.logP)
            except:
                pass

        # Get molecule image
        if hasattr(m, 'image') and m.image is not None:
            try:
                buf = io.BytesIO()
                m.image.save(buf, format='PNG')
                info.image_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')
            except Exception as e:
                pass

        info.validated = True

    except Exception as e:
        info.error = str(e)
        info.validated = False

    return info


def compute_triangular_pdf(c_min: float, c_mode: float, c_max: float, n_points: int = 200) -> Dict[str, List[float]]:
    """
    Compute triangular PDF for visualization.

    Parameters
    ----------
    c_min, c_mode, c_max : float
        Triangular distribution parameters.
    n_points : int
        Number of points for plotting.

    Returns
    -------
    dict
        {'x': [...], 'y': [...]} for plotting.
    """
    if c_max <= c_min:
        c_max = c_min + 1.0
    if c_mode < c_min:
        c_mode = c_min
    if c_mode > c_max:
        c_mode = c_max

    x = np.linspace(c_min, c_max, n_points)
    y = np.zeros_like(x)

    # Height at mode
    h = 2.0 / (c_max - c_min)

    # Left side: linear rise from c_min to c_mode
    mask_left = (x >= c_min) & (x <= c_mode)
    if c_mode > c_min:
        y[mask_left] = h * (x[mask_left] - c_min) / (c_mode - c_min)

    # Right side: linear fall from c_mode to c_max
    mask_right = (x > c_mode) & (x <= c_max)
    if c_max > c_mode:
        y[mask_right] = h * (c_max - x[mask_right]) / (c_max - c_mode)

    return {'x': x.tolist(), 'y': y.tolist()}


# =============================================================================
# FastAPI Application
# =============================================================================

app = FastAPI(
    title="SFPPy Family Editor",
    description="Interactive editor for substance families",
    version="1.0.0",
)

# Mount static files
STATIC_DIR = Path(__file__).parent / "static"
TEMPLATES_DIR = Path(__file__).parent / "templates"

app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")
templates = Jinja2Templates(directory=str(TEMPLATES_DIR))


# =============================================================================
# Routes
# =============================================================================

@app.get("/", response_class=HTMLResponse)
async def home(request: Request):
    """Main page - family editor."""
    return templates.TemplateResponse("index.html", {
        "request": request,
        "family": current_family,
    })


@app.post("/api/lookup")
async def api_lookup(identifier: str = Form(...)):
    """Look up substance by name or CAS."""
    info = lookup_substance(identifier)
    return JSONResponse({
        "success": info.validated,
        "substance": asdict(info),
    })


@app.post("/api/add-substance")
async def api_add_substance(
    identifier: str = Form(...),
    c0_min: float = Form(0.0),
    c0_likely: float = Form(100.0),
    c0_max: float = Form(500.0),
    notes: str = Form(""),
):
    """Add substance to current family."""
    global current_family

    info = lookup_substance(identifier)
    info.C0_min = c0_min
    info.C0_likely = c0_likely
    info.C0_max = c0_max
    info.notes = notes

    # Check for duplicates
    for existing in current_family.substances:
        if existing.cas and info.cas and existing.cas == info.cas:
            return JSONResponse({
                "success": False,
                "error": f"Substance with CAS {info.cas} already exists in family",
            })
        if existing.name.lower() == info.name.lower():
            return JSONResponse({
                "success": False,
                "error": f"Substance '{info.name}' already exists in family",
            })

    current_family.substances.append(info)

    return JSONResponse({
        "success": True,
        "substance": asdict(info),
        "family_size": len(current_family.substances),
    })


@app.post("/api/remove-substance")
async def api_remove_substance(index: int = Form(...)):
    """Remove substance from current family by index."""
    global current_family

    if 0 <= index < len(current_family.substances):
        removed = current_family.substances.pop(index)
        return JSONResponse({
            "success": True,
            "removed": removed.name,
            "family_size": len(current_family.substances),
        })

    return JSONResponse({
        "success": False,
        "error": "Invalid index",
    })


@app.post("/api/update-substance")
async def api_update_substance(
    index: int = Form(...),
    c0_min: float = Form(0.0),
    c0_likely: float = Form(100.0),
    c0_max: float = Form(500.0),
    notes: str = Form(""),
):
    """Update substance concentration parameters."""
    global current_family

    if 0 <= index < len(current_family.substances):
        current_family.substances[index].C0_min = c0_min
        current_family.substances[index].C0_likely = c0_likely
        current_family.substances[index].C0_max = c0_max
        current_family.substances[index].notes = notes
        return JSONResponse({"success": True})

    return JSONResponse({"success": False, "error": "Invalid index"})


@app.post("/api/update-family")
async def api_update_family(
    name: str = Form(...),
    description: str = Form(""),
):
    """Update family metadata."""
    global current_family
    current_family.name = name
    current_family.description = description
    return JSONResponse({"success": True})


@app.get("/api/family")
async def api_get_family():
    """Get current family data."""
    return JSONResponse({
        "name": current_family.name,
        "description": current_family.description,
        "substances": [asdict(s) for s in current_family.substances],
    })


@app.post("/api/new-family")
async def api_new_family():
    """Create new empty family."""
    global current_family
    current_family = FamilyInfo()
    return JSONResponse({"success": True})


@app.get("/api/export-yaml")
async def api_export_yaml():
    """Export family as YAML."""
    yaml_content = current_family.to_yaml()
    return JSONResponse({
        "yaml": yaml_content,
        "filename": f"{current_family.name}.yml",
    })


@app.post("/api/import-yaml")
async def api_import_yaml(file: UploadFile = File(...)):
    """Import family from YAML file."""
    global current_family

    try:
        content = await file.read()
        data = yaml.safe_load(content.decode('utf-8'))

        current_family = FamilyInfo(
            name=data.get('name', 'imported_family'),
            description=data.get('description', ''),
        )

        # Process substances
        errors = []
        for s_data in data.get('substances', []):
            identifier = s_data.get('cas') or s_data.get('name')
            if not identifier:
                continue

            info = lookup_substance(identifier)
            info.C0_min = float(s_data.get('C0_min', 0.0))
            info.C0_likely = float(s_data.get('C0_likely', 100.0))
            info.C0_max = float(s_data.get('C0_max', 500.0))
            info.notes = s_data.get('notes', '')

            if info.error:
                errors.append(f"{identifier}: {info.error}")

            current_family.substances.append(info)

        return JSONResponse({
            "success": True,
            "family_name": current_family.name,
            "substances_count": len(current_family.substances),
            "errors": errors,
        })

    except Exception as e:
        return JSONResponse({
            "success": False,
            "error": str(e),
        })


@app.get("/api/distribution/{index}")
async def api_get_distribution(index: int):
    """Get triangular distribution data for plotting."""
    if 0 <= index < len(current_family.substances):
        s = current_family.substances[index]
        dist = compute_triangular_pdf(s.C0_min, s.C0_likely, s.C0_max)
        return JSONResponse({
            "success": True,
            "distribution": dist,
            "params": {
                "min": s.C0_min,
                "likely": s.C0_likely,
                "max": s.C0_max,
            }
        })
    return JSONResponse({"success": False, "error": "Invalid index"})


@app.get("/api/all-distributions")
async def api_all_distributions():
    """Get all substance distributions for comparison chart."""
    distributions = []
    for i, s in enumerate(current_family.substances):
        dist = compute_triangular_pdf(s.C0_min, s.C0_likely, s.C0_max)
        distributions.append({
            "name": s.name,
            "index": i,
            "validated": s.validated,
            "distribution": dist,
            "params": {
                "min": s.C0_min,
                "likely": s.C0_likely,
                "max": s.C0_max,
            }
        })
    return JSONResponse({"distributions": distributions})


@app.get("/api/validate-all")
async def api_validate_all():
    """Re-validate all substances in current family."""
    global current_family

    results = []
    for i, s in enumerate(current_family.substances):
        if not s.validated:
            # Try to re-validate
            identifier = s.cas or s.name
            new_info = lookup_substance(identifier)
            if new_info.validated:
                # Preserve C0 values
                new_info.C0_min = s.C0_min
                new_info.C0_likely = s.C0_likely
                new_info.C0_max = s.C0_max
                new_info.notes = s.notes
                current_family.substances[i] = new_info
                results.append({"index": i, "name": s.name, "status": "fixed"})
            else:
                results.append({"index": i, "name": s.name, "status": "failed", "error": new_info.error})
        else:
            results.append({"index": i, "name": s.name, "status": "ok"})

    return JSONResponse({"results": results})


# =============================================================================
# Download YAML file
# =============================================================================

@app.get("/download-yaml")
async def download_yaml():
    """Download family as YAML file."""
    yaml_content = current_family.to_yaml()

    # Create temporary file
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as f:
        f.write(yaml_content)
        temp_path = f.name

    return FileResponse(
        temp_path,
        media_type='application/x-yaml',
        filename=f"{current_family.name}.yml",
    )


# =============================================================================
# Working Folder & Discovery
# =============================================================================

# Default working folder for families
WORKING_FOLDER: Path = Path(__file__).parent.parent / "examples" / "families"


@app.get("/api/working-folder")
async def api_get_working_folder():
    """Get current working folder."""
    return JSONResponse({
        "path": str(WORKING_FOLDER),
        "exists": WORKING_FOLDER.exists(),
    })


@app.post("/api/set-working-folder")
async def api_set_working_folder(path: str = Form(...)):
    """Set working folder for family discovery."""
    global WORKING_FOLDER
    new_path = Path(path).expanduser().resolve()

    if not new_path.exists():
        return JSONResponse({
            "success": False,
            "error": f"Path does not exist: {new_path}",
        })

    if not new_path.is_dir():
        return JSONResponse({
            "success": False,
            "error": f"Path is not a directory: {new_path}",
        })

    WORKING_FOLDER = new_path
    return JSONResponse({
        "success": True,
        "path": str(WORKING_FOLDER),
    })


@app.get("/api/discover-families")
async def api_discover_families():
    """Discover family YAML files in working folder."""
    if not WORKING_FOLDER.exists():
        return JSONResponse({
            "success": False,
            "error": f"Working folder does not exist: {WORKING_FOLDER}",
        })

    families = []
    for yml_file in WORKING_FOLDER.glob("*.yml"):
        try:
            with open(yml_file, 'r', encoding='utf-8') as f:
                data = yaml.safe_load(f)

            # Check if it's a valid family file
            if not isinstance(data, dict):
                continue
            if 'substances' not in data:
                continue

            families.append({
                "filename": yml_file.name,
                "path": str(yml_file),
                "name": data.get('name', yml_file.stem),
                "description": data.get('description', ''),
                "substance_count": len(data.get('substances', [])),
            })
        except Exception as e:
            families.append({
                "filename": yml_file.name,
                "path": str(yml_file),
                "error": str(e),
            })

    return JSONResponse({
        "success": True,
        "folder": str(WORKING_FOLDER),
        "families": families,
    })


@app.post("/api/load-family-from-folder")
async def api_load_family_from_folder(filename: str = Form(...)):
    """Load a family from the working folder."""
    global current_family

    file_path = WORKING_FOLDER / filename
    if not file_path.exists():
        return JSONResponse({
            "success": False,
            "error": f"File not found: {filename}",
        })

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            data = yaml.safe_load(f)

        current_family = FamilyInfo(
            name=data.get('name', file_path.stem),
            description=data.get('description', ''),
        )

        # Process substances with validation
        errors = []
        validated_count = 0
        for s_data in data.get('substances', []):
            identifier = s_data.get('cas') or s_data.get('name')
            if not identifier:
                continue

            info = lookup_substance(identifier)
            info.C0_min = float(s_data.get('C0_min', 0.0))
            info.C0_likely = float(s_data.get('C0_likely', 100.0))
            info.C0_max = float(s_data.get('C0_max', 500.0))
            info.notes = s_data.get('notes', '')

            if info.error:
                errors.append(f"{identifier}: {info.error}")
            else:
                validated_count += 1

            current_family.substances.append(info)

        return JSONResponse({
            "success": True,
            "family_name": current_family.name,
            "substances_count": len(current_family.substances),
            "validated_count": validated_count,
            "errors": errors,
        })

    except Exception as e:
        return JSONResponse({
            "success": False,
            "error": str(e),
        })


@app.post("/api/save-family-to-folder")
async def api_save_family_to_folder(filename: str = Form(None)):
    """Save current family to working folder."""
    if not WORKING_FOLDER.exists():
        WORKING_FOLDER.mkdir(parents=True, exist_ok=True)

    if filename is None:
        filename = f"{current_family.name}.yml"

    if not filename.endswith('.yml'):
        filename += '.yml'

    file_path = WORKING_FOLDER / filename

    try:
        yaml_content = current_family.to_yaml()
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(yaml_content)

        return JSONResponse({
            "success": True,
            "path": str(file_path),
            "filename": filename,
        })

    except Exception as e:
        return JSONResponse({
            "success": False,
            "error": str(e),
        })


@app.get("/api/validate-family-file")
async def api_validate_family_file(filename: str):
    """Validate a family file from the working folder without loading it."""
    file_path = WORKING_FOLDER / filename
    if not file_path.exists():
        return JSONResponse({
            "success": False,
            "error": f"File not found: {filename}",
        })

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            data = yaml.safe_load(f)

        if not isinstance(data, dict):
            return JSONResponse({
                "success": False,
                "error": "Invalid YAML structure",
            })

        if 'substances' not in data:
            return JSONResponse({
                "success": False,
                "error": "Missing 'substances' key",
            })

        # Validate each substance
        results = []
        for s_data in data.get('substances', []):
            identifier = s_data.get('cas') or s_data.get('name')
            if not identifier:
                results.append({
                    "identifier": "unknown",
                    "status": "error",
                    "error": "Missing name or CAS",
                })
                continue

            info = lookup_substance(identifier)
            if info.validated:
                results.append({
                    "identifier": identifier,
                    "name": info.name,
                    "cas": info.cas,
                    "mass_g_mol": info.mass_g_mol,
                    "status": "valid",
                })
            else:
                results.append({
                    "identifier": identifier,
                    "status": "error",
                    "error": info.error,
                })

        valid_count = sum(1 for r in results if r.get('status') == 'valid')
        error_count = sum(1 for r in results if r.get('status') == 'error')

        return JSONResponse({
            "success": True,
            "filename": filename,
            "family_name": data.get('name', filename),
            "total": len(results),
            "valid": valid_count,
            "errors": error_count,
            "substances": results,
        })

    except Exception as e:
        return JSONResponse({
            "success": False,
            "error": str(e),
        })


# =============================================================================
# Health check
# =============================================================================

@app.get("/health")
async def health():
    """Health check endpoint."""
    return {
        "status": "ok",
        "family_substances": len(current_family.substances),
        "working_folder": str(WORKING_FOLDER),
    }


# =============================================================================
# Run with uvicorn
# =============================================================================

if __name__ == "__main__":
    import uvicorn
    print(f"""
╔══════════════════════════════════════════════════════════════╗
║           SFPPy Family Editor — Starting Server              ║
╚══════════════════════════════════════════════════════════════╝

Open your browser at: http://localhost:8000
""")
    uvicorn.run(app, host="0.0.0.0", port=8000)
