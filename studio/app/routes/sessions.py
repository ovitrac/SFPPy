"""
Session Routes - Simulation Session Management

API endpoints for managing simulation sessions with substance storage.
"""

import sys
from pathlib import Path
from typing import Optional, List

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

# Add parent paths
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from studio.app.session import (
    get_session_store,
    load_migrant_toxtree,
    create_unknown_substance,
    compute_substance_properties,
    convert_concentration,
    get_concentration_for_computation,
    CONCENTRATION_UNITS,
    SubstanceData,
    LayerSubstanceAssignment,
    ContactStep,
    LayerConfig,
    TTC_VALUES,
    CFTTC_VALUES,
)

router = APIRouter()


# ========== DATA MODELS ==========

class CreateSessionInput(BaseModel):
    name: str = Field(default="Untitled Simulation")


class AddSubstanceInput(BaseModel):
    query: str = Field(..., description="Substance name, CAS, or CID")


class CreateUnknownSubstanceInput(BaseModel):
    name: str = Field(default="unknown")
    mw: float = Field(default=500.0, gt=0, description="Molecular weight (g/mol)")
    logP: float = Field(default=5.0, description="Log partition coefficient")


class AssignSubstanceInput(BaseModel):
    substance_id: str = Field(..., description="Substance ID in session")
    layer_index: int = Field(..., ge=1, le=10, description="Layer index")
    C0: float = Field(..., gt=0, description="Initial concentration")
    C0_unit: str = Field(default="mg/kg", description="Concentration unit")
    D_override: Optional[float] = Field(None, description="Override diffusivity (m²/s)")
    k_override: Optional[float] = Field(None, description="Override partition coefficient")


class UpdateAssignmentInput(BaseModel):
    C0: Optional[float] = Field(None, gt=0)
    C0_unit: Optional[str] = None
    D_override: Optional[float] = None
    k_override: Optional[float] = None


class AddContactStepInput(BaseModel):
    temperature: float = Field(..., description="Temperature value")
    temperature_unit: str = Field(default="degC", description="Temperature unit")
    duration: float = Field(..., gt=0, description="Duration value")
    duration_unit: str = Field(default="days", description="Duration unit")
    with_food: bool = Field(default=True, description="Contact with food (False for set-off)")
    simulant: Optional[str] = Field(None, description="Food simulant code")


class AddLayerInput(BaseModel):
    polymer: str = Field(..., description="Polymer code")
    thickness: float = Field(..., gt=0, description="Thickness value")
    thickness_unit: str = Field(default="um", description="Thickness unit")


# ========== SESSION ENDPOINTS ==========

@router.post("/create")
async def create_session(data: CreateSessionInput):
    """Create a new simulation session."""
    store = get_session_store()
    session = store.create_session(data.name)

    return JSONResponse({
        "success": True,
        "session": {
            "id": session.id,
            "name": session.name,
            "created_at": session.created_at.isoformat(),
        },
        "message": f"Session '{session.name}' created",
    })


@router.get("/list")
async def list_sessions():
    """List all active sessions."""
    store = get_session_store()
    sessions = store.list_sessions()

    return JSONResponse({
        "success": True,
        "sessions": sessions,
        "count": len(sessions),
    })


@router.get("/{session_id}")
async def get_session(session_id: str):
    """Get full session state."""
    store = get_session_store()
    session = store.get_session(session_id)

    if session is None:
        raise HTTPException(status_code=404, detail="Session not found or expired")

    return JSONResponse({
        "success": True,
        "session": session.to_dict(),
    })


@router.delete("/{session_id}")
async def delete_session(session_id: str):
    """Delete a session."""
    store = get_session_store()

    if store.delete_session(session_id):
        return JSONResponse({
            "success": True,
            "message": "Session deleted",
        })
    else:
        raise HTTPException(status_code=404, detail="Session not found")


@router.get("/{session_id}/validate")
async def validate_session(session_id: str):
    """Validate session is ready for simulation."""
    store = get_session_store()
    session = store.get_session(session_id)

    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    is_valid, errors, warnings = session.is_valid_for_simulation()

    return JSONResponse({
        "success": True,
        "valid": is_valid,
        "errors": errors,
        "warnings": warnings,
        "summary": {
            "layers_count": len(session.layers),
            "substances_count": len(session.substances),
            "assignments_count": len(session.assignments),
            "steps_count": len(session.contact_steps),
        },
    })


# ========== SUBSTANCE ENDPOINTS ==========

@router.get("/{session_id}/substances")
async def list_session_substances(session_id: str):
    """List all substances in session."""
    store = get_session_store()
    session = store.get_session(session_id)

    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    substances = [s.to_dict() for s in session.substances.values()]

    return JSONResponse({
        "success": True,
        "substances": substances,
        "count": len(substances),
    })


@router.post("/{session_id}/substances")
async def add_substance_to_session(session_id: str, data: AddSubstanceInput):
    """
    Add a substance to the session.

    Loads substance data via migrantToxtree, including:
    - PubChem properties (MW, logP, SMILES, etc.)
    - CAS array (multiple CAS numbers if available)
    - Cramer classification
    - TTC and CF_TTC values
    - Regulatory status (EU, US, CN)
    """
    store = get_session_store()
    session = store.get_session(session_id)

    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    # Load substance via migrantToxtree
    substance = load_migrant_toxtree(data.query)

    if substance is None:
        raise HTTPException(
            status_code=404,
            detail=f"Substance '{data.query}' not found in PubChem"
        )

    # Check if already in session
    if substance.id in session.substances:
        return JSONResponse({
            "success": True,
            "substance": substance.to_dict(),
            "message": "Substance already in session",
            "added": False,
        })

    # Add to session
    session.substances[substance.id] = substance
    store.update_session(session)

    return JSONResponse({
        "success": True,
        "substance": substance.to_dict(),
        "message": f"Substance '{substance.name}' added to session",
        "added": True,
    })


@router.post("/{session_id}/substances/unknown")
async def add_unknown_substance(session_id: str, data: CreateUnknownSubstanceInput):
    """
    Add an 'unknown' substance with default properties.

    Used when no specific substance data is available.
    Default: D=1e-12 m²/s, k=1
    """
    store = get_session_store()
    session = store.get_session(session_id)

    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    substance = create_unknown_substance(data.name, data.mw, data.logP)

    # Check if already in session
    if substance.id in session.substances:
        return JSONResponse({
            "success": True,
            "substance": substance.to_dict(),
            "message": "Unknown substance already in session",
            "added": False,
        })

    session.substances[substance.id] = substance
    store.update_session(session)

    return JSONResponse({
        "success": True,
        "substance": substance.to_dict(),
        "message": f"Unknown substance '{substance.name}' added with D=1e-12, k=1",
        "added": True,
    })


@router.delete("/{session_id}/substances/{substance_id}")
async def remove_substance_from_session(session_id: str, substance_id: str):
    """Remove a substance from the session."""
    store = get_session_store()
    session = store.get_session(session_id)

    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    if substance_id not in session.substances:
        raise HTTPException(status_code=404, detail="Substance not in session")

    # Remove substance
    del session.substances[substance_id]

    # Remove any assignments for this substance
    session.assignments = [
        a for a in session.assignments
        if a.substance_id != substance_id
    ]

    store.update_session(session)

    return JSONResponse({
        "success": True,
        "message": f"Substance '{substance_id}' removed from session",
    })


@router.get("/{session_id}/substances/{substance_id}")
async def get_substance_details(session_id: str, substance_id: str):
    """Get full details for a substance in session."""
    store = get_session_store()
    session = store.get_session(session_id)

    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    if substance_id not in session.substances:
        raise HTTPException(status_code=404, detail="Substance not in session")

    substance = session.substances[substance_id]

    return JSONResponse({
        "success": True,
        "substance": substance.to_dict(),
    })


@router.get("/{session_id}/substances/{substance_id}/properties")
async def compute_substance_props(
    session_id: str,
    substance_id: str,
    polymer: str = Query(..., description="Polymer code"),
    temperature_C: float = Query(25.0, description="Temperature in Celsius"),
):
    """
    Compute D and k for a substance in a specific polymer.

    Uses SFPPy layer module for proper computation.
    """
    store = get_session_store()
    session = store.get_session(session_id)

    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    if substance_id not in session.substances:
        raise HTTPException(status_code=404, detail="Substance not in session")

    substance = session.substances[substance_id]
    props = compute_substance_properties(substance, polymer, temperature_C)

    return JSONResponse({
        "success": True,
        "substance_id": substance_id,
        "polymer": polymer,
        "temperature_C": temperature_C,
        "properties": props,
    })


# ========== ASSIGNMENT ENDPOINTS ==========

@router.post("/{session_id}/assignments")
async def assign_substance_to_layer(session_id: str, data: AssignSubstanceInput):
    """
    Assign a substance to a layer with initial concentration.

    Automatically computes D and k based on layer polymer and temperature.
    """
    store = get_session_store()
    session = store.get_session(session_id)

    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    if data.substance_id not in session.substances:
        raise HTTPException(status_code=404, detail="Substance not in session")

    # Check if layer exists
    layer_exists = any(l.index == data.layer_index for l in session.layers)
    if not layer_exists and len(session.layers) > 0:
        raise HTTPException(status_code=400, detail=f"Layer {data.layer_index} not found")

    # Check for existing assignment
    existing = next(
        (a for a in session.assignments
         if a.substance_id == data.substance_id and a.layer_index == data.layer_index),
        None
    )

    if existing:
        # Update existing assignment
        existing.C0 = data.C0
        existing.C0_unit = data.C0_unit
        existing.D_override = data.D_override
        existing.k_override = data.k_override
        message = "Assignment updated"
    else:
        # Create new assignment
        assignment = LayerSubstanceAssignment(
            substance_id=data.substance_id,
            layer_index=data.layer_index,
            C0=data.C0,
            C0_unit=data.C0_unit,
            D_override=data.D_override,
            k_override=data.k_override,
        )
        session.assignments.append(assignment)
        existing = assignment
        message = "Substance assigned to layer"

    # Compute properties if layer is defined
    layer = next((l for l in session.layers if l.index == data.layer_index), None)
    if layer:
        substance = session.substances[data.substance_id]
        props = compute_substance_properties(substance, layer.polymer, 25.0)  # Default T

        existing.D_computed = props["D"]
        existing.k_computed = props["k"]
        existing.D_final = data.D_override if data.D_override else props["D"]
        existing.k_final = data.k_override if data.k_override else props["k"]

    store.update_session(session)

    return JSONResponse({
        "success": True,
        "assignment": existing.to_dict(),
        "message": message,
    })


@router.get("/{session_id}/assignments")
async def list_assignments(session_id: str):
    """List all substance-to-layer assignments."""
    store = get_session_store()
    session = store.get_session(session_id)

    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    assignments = [a.to_dict() for a in session.assignments]

    return JSONResponse({
        "success": True,
        "assignments": assignments,
        "count": len(assignments),
    })


@router.delete("/{session_id}/assignments/{substance_id}/{layer_index}")
async def remove_assignment(session_id: str, substance_id: str, layer_index: int):
    """Remove a substance assignment from a layer."""
    store = get_session_store()
    session = store.get_session(session_id)

    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    original_count = len(session.assignments)
    session.assignments = [
        a for a in session.assignments
        if not (a.substance_id == substance_id and a.layer_index == layer_index)
    ]

    if len(session.assignments) == original_count:
        raise HTTPException(status_code=404, detail="Assignment not found")

    store.update_session(session)

    return JSONResponse({
        "success": True,
        "message": "Assignment removed",
    })


# ========== LAYER ENDPOINTS ==========

@router.post("/{session_id}/layers")
async def add_layer(session_id: str, data: AddLayerInput):
    """Add a layer to the session."""
    store = get_session_store()
    session = store.get_session(session_id)

    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    # Auto-assign index
    new_index = len(session.layers) + 1

    layer = LayerConfig(
        index=new_index,
        polymer=data.polymer,
        thickness=data.thickness,
        thickness_unit=data.thickness_unit,
    )

    session.layers.append(layer)
    store.update_session(session)

    return JSONResponse({
        "success": True,
        "layer": layer.to_dict(),
        "message": f"Layer {new_index} added",
    })


@router.get("/{session_id}/layers")
async def list_layers(session_id: str):
    """List all layers in session."""
    store = get_session_store()
    session = store.get_session(session_id)

    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    layers = [l.to_dict() for l in session.layers]

    return JSONResponse({
        "success": True,
        "layers": layers,
        "count": len(layers),
    })


@router.delete("/{session_id}/layers/{layer_index}")
async def remove_layer(session_id: str, layer_index: int):
    """Remove a layer from session."""
    store = get_session_store()
    session = store.get_session(session_id)

    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    original_count = len(session.layers)
    session.layers = [l for l in session.layers if l.index != layer_index]

    if len(session.layers) == original_count:
        raise HTTPException(status_code=404, detail="Layer not found")

    # Re-index remaining layers
    for i, layer in enumerate(session.layers, start=1):
        layer.index = i

    # Remove assignments for deleted layer
    session.assignments = [a for a in session.assignments if a.layer_index != layer_index]

    store.update_session(session)

    return JSONResponse({
        "success": True,
        "message": f"Layer {layer_index} removed",
    })


# ========== CONTACT STEP ENDPOINTS ==========

@router.post("/{session_id}/steps")
async def add_contact_step(session_id: str, data: AddContactStepInput):
    """Add a contact step to the session."""
    store = get_session_store()
    session = store.get_session(session_id)

    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    # Auto-assign index
    new_index = len(session.contact_steps) + 1

    step = ContactStep(
        index=new_index,
        temperature=data.temperature,
        temperature_unit=data.temperature_unit,
        duration=data.duration,
        duration_unit=data.duration_unit,
        with_food=data.with_food,
        simulant=data.simulant,
    )

    session.contact_steps.append(step)
    store.update_session(session)

    return JSONResponse({
        "success": True,
        "step": step.to_dict(),
        "message": f"Contact step {new_index} added",
    })


@router.get("/{session_id}/steps")
async def list_contact_steps(session_id: str):
    """List all contact steps in session."""
    store = get_session_store()
    session = store.get_session(session_id)

    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    steps = [s.to_dict() for s in session.contact_steps]

    return JSONResponse({
        "success": True,
        "steps": steps,
        "count": len(steps),
    })


@router.delete("/{session_id}/steps/{step_index}")
async def remove_contact_step(session_id: str, step_index: int):
    """Remove a contact step from session."""
    store = get_session_store()
    session = store.get_session(session_id)

    if session is None:
        raise HTTPException(status_code=404, detail="Session not found")

    original_count = len(session.contact_steps)
    session.contact_steps = [s for s in session.contact_steps if s.index != step_index]

    if len(session.contact_steps) == original_count:
        raise HTTPException(status_code=404, detail="Step not found")

    # Re-index remaining steps
    for i, step in enumerate(session.contact_steps, start=1):
        step.index = i

    store.update_session(session)

    return JSONResponse({
        "success": True,
        "message": f"Contact step {step_index} removed",
    })


# ========== CONCENTRATION UNITS ==========

@router.get("/units/concentration")
async def list_concentration_units():
    """List available concentration units."""
    return JSONResponse({
        "success": True,
        "units": CONCENTRATION_UNITS,
        "default": "mg/kg",
    })


@router.post("/units/convert")
async def convert_concentration_endpoint(
    value: float = Query(..., description="Value to convert"),
    from_unit: str = Query(..., description="Source unit"),
    to_unit: str = Query(..., description="Target unit"),
    polymer_density: float = Query(1000.0, description="Polymer density (kg/m³)"),
):
    """Convert concentration between units."""
    try:
        result = convert_concentration(value, from_unit, to_unit, polymer_density)
        return JSONResponse({
            "success": True,
            "input": {"value": value, "unit": from_unit},
            "output": {"value": result, "unit": to_unit},
            "polymer_density": polymer_density,
        })
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


# ========== TTC REFERENCE ==========

@router.get("/reference/ttc")
async def get_ttc_reference():
    """Get TTC (Threshold of Toxicological Concern) reference values."""
    return JSONResponse({
        "success": True,
        "ttc": {
            "values": TTC_VALUES,
            "unit": "µg/kg bw/day",
            "description": "Threshold of Toxicological Concern by Cramer class",
            "classes": ["Unknown (0)", "Class I", "Class II", "Class III"],
        },
        "cf_ttc": {
            "values": CFTTC_VALUES,
            "unit": "mg/kg food",
            "description": "Conversion factor for food intake (60 kg bw, 1 kg food/day)",
        },
    })


# ========== DIRECT D/k/k0 COMPUTATION ==========

@router.get("/compute/diffusivity")
async def compute_diffusivity_direct(
    substance: str = Query(..., description="Substance name, CAS, or CID"),
    polymer: str = Query(..., description="Polymer code (e.g., LDPE, PP, PET)"),
    temperature_C: float = Query(40.0, description="Temperature in Celsius"),
):
    """
    Compute D, k, D0, k0 for a substance in a specific polymer.

    This is a direct computation endpoint that doesn't require a session.
    Uses SFPPy layer module for proper Piringer/FHP calculations.

    Returns:
    - D: Diffusivity at specified temperature (m²/s)
    - D0: Reference diffusivity (m²/s)
    - k: Partition coefficient at specified temperature
    - k0: Reference partition coefficient
    - method: Computation method used ('sfppy' or 'estimation')
    """
    # Load substance - try multiple name variants
    substance_data = load_migrant_toxtree(substance)

    # If not found and name contains parentheses, try the short name
    if substance_data is None and '(' in substance:
        short_name = substance.split('(')[0].strip()
        substance_data = load_migrant_toxtree(short_name)

    # If still not found, try common aliases for known substances
    if substance_data is None:
        SUBSTANCE_ALIASES = {
            'bht': '128-37-0',  # CAS for BHT
            'butylated hydroxytoluene': '128-37-0',
            'irganox 1076': '2082-79-3',
            'irganox 1010': '6683-19-8',
            'deha': '103-23-1',
            'dehp': '117-81-7',
        }
        alias_query = SUBSTANCE_ALIASES.get(substance.lower().split('(')[0].strip())
        if alias_query:
            substance_data = load_migrant_toxtree(alias_query)

    if substance_data is None:
        raise HTTPException(
            status_code=404,
            detail=f"Substance '{substance}' not found in PubChem. Use a valid CAS number, IUPAC name, brand name (e.g., Irganox 1076), or known acronym (e.g., BHT)."
        )

    # Compute properties
    props = compute_substance_properties(substance_data, polymer, temperature_C)

    return JSONResponse({
        "success": True,
        "substance": {
            "name": substance_data.name,
            "cid": substance_data.cid,
            "mw": substance_data.mw,
            "logP": substance_data.logP,
        },
        "polymer": polymer,
        "temperature_C": temperature_C,
        "properties": {
            "D": props["D"],
            "D_unit": "m²/s",
            "D0": props["D0"],
            "D0_unit": "m²/s",
            "k": props["k"],
            "k_unit": "dimensionless",
            "k0": props["k0"],
            "k0_unit": "dimensionless",
        },
        "method": props["method"],
    })


@router.get("/compute/diffusivity/batch")
async def compute_diffusivity_batch(
    substance: str = Query(..., description="Substance name, CAS, or CID"),
    polymers: str = Query("LDPE,PP,PET", description="Comma-separated polymer codes"),
    temperature_C: float = Query(40.0, description="Temperature in Celsius"),
):
    """
    Compute D and k for a substance across multiple polymers.

    Useful for comparing migration behavior in different materials.
    """
    # Load substance
    substance_data = load_migrant_toxtree(substance)

    if substance_data is None:
        raise HTTPException(
            status_code=404,
            detail=f"Substance '{substance}' not found in PubChem. Use a valid CAS number, IUPAC name, brand name (e.g., Irganox 1076), or known acronym (e.g., BHT)."
        )

    # Parse polymer list
    polymer_list = [p.strip() for p in polymers.split(",")]

    results = []
    for polymer in polymer_list:
        props = compute_substance_properties(substance_data, polymer, temperature_C)
        results.append({
            "polymer": polymer,
            "D": props["D"],
            "D0": props["D0"],
            "k": props["k"],
            "k0": props["k0"],
            "method": props["method"],
        })

    return JSONResponse({
        "success": True,
        "substance": {
            "name": substance_data.name,
            "cid": substance_data.cid,
            "mw": substance_data.mw,
            "logP": substance_data.logP,
        },
        "temperature_C": temperature_C,
        "results": results,
    })


@router.get("/compute/diffusivity/temperature")
async def compute_diffusivity_temperature_range(
    substance: str = Query(..., description="Substance name, CAS, or CID"),
    polymer: str = Query("LDPE", description="Polymer code"),
    temp_min: float = Query(4.0, description="Minimum temperature (°C)"),
    temp_max: float = Query(60.0, description="Maximum temperature (°C)"),
    temp_step: float = Query(10.0, description="Temperature step (°C)"),
):
    """
    Compute D at multiple temperatures for Arrhenius analysis.

    Useful for understanding temperature dependence of diffusion.
    """
    # Load substance
    substance_data = load_migrant_toxtree(substance)

    if substance_data is None:
        raise HTTPException(
            status_code=404,
            detail=f"Substance '{substance}' not found in PubChem. Use a valid CAS number, IUPAC name, brand name (e.g., Irganox 1076), or known acronym (e.g., BHT)."
        )

    # Generate temperature range
    temperatures = []
    t = temp_min
    while t <= temp_max:
        temperatures.append(t)
        t += temp_step

    results = []
    for T in temperatures:
        props = compute_substance_properties(substance_data, polymer, T)
        results.append({
            "temperature_C": T,
            "temperature_K": T + 273.15,
            "D": props["D"],
            "k": props["k"],
        })

    return JSONResponse({
        "success": True,
        "substance": {
            "name": substance_data.name,
            "mw": substance_data.mw,
        },
        "polymer": polymer,
        "results": results,
    })


# =============================================================================
# SESSION FILE ENDPOINTS (NEW - .sfppy.json format)
# =============================================================================

from studio.app.models.session import Session as SessionFile
from studio.app.utils.session_io import (
    load_session_file,
    save_session_file,
    list_session_files,
    validate_session,
    export_session_summary,
    session_to_patankar_inputs,
)


class LoadSessionFileInput(BaseModel):
    filepath: str = Field(..., description="Path to .sfppy.json file")


class SaveSessionFileInput(BaseModel):
    filepath: str = Field(..., description="Target file path")
    overwrite: bool = Field(default=False, description="Overwrite if exists")


@router.get("/files/list")
async def list_example_sessions(
    directory: str = Query(None, description="Directory to list (default: examples)")
):
    """
    List available session files (.sfppy.json).

    Returns a list of session files with metadata.
    """
    # Default to examples directory
    if directory is None:
        directory = str(Path(__file__).parent.parent.parent / "examples")

    files = list_session_files(directory)

    return JSONResponse({
        "success": True,
        "directory": directory,
        "files": files,
        "count": len(files),
    })


@router.get("/files/examples")
async def list_builtin_examples():
    """
    List built-in example session files.

    These are the reference examples for testing and demonstration.
    """
    examples_dir = Path(__file__).parent.parent.parent / "examples"
    files = list_session_files(str(examples_dir))

    return JSONResponse({
        "success": True,
        "examples_directory": str(examples_dir),
        "examples": files,
        "count": len(files),
    })


@router.post("/files/load")
async def load_session_from_file(data: LoadSessionFileInput):
    """
    Load a session from a .sfppy.json file.

    Returns the full session data including substances, layers, steps.
    """
    try:
        session = load_session_file(data.filepath)
        is_valid, errors, warnings = validate_session(session)

        return JSONResponse({
            "success": True,
            "filepath": data.filepath,
            "session": session.model_dump(mode='json', exclude_none=True),
            "validation": {
                "valid": is_valid,
                "errors": errors,
                "warnings": warnings,
            },
            "summary": export_session_summary(session),
        })
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Failed to load session: {e}")


@router.get("/files/load/{example_name}")
async def load_example_session(example_name: str):
    """
    Load a built-in example session by name.

    Example names: example1, example1_extensions, example3, example3_variant
    """
    examples_dir = Path(__file__).parent.parent.parent / "examples"

    # Try with and without .sfppy.json extension
    filepath = examples_dir / f"{example_name}.sfppy.json"
    if not filepath.exists():
        filepath = examples_dir / example_name
        if not filepath.exists():
            raise HTTPException(
                status_code=404,
                detail=f"Example '{example_name}' not found. Available: example1, example1_extensions, example3, example3_variant"
            )

    try:
        session = load_session_file(str(filepath))
        is_valid, errors, warnings = validate_session(session)

        return JSONResponse({
            "success": True,
            "example_name": example_name,
            "filepath": str(filepath),
            "session": session.model_dump(mode='json', exclude_none=True),
            "validation": {
                "valid": is_valid,
                "errors": errors,
                "warnings": warnings,
            },
            "summary": export_session_summary(session),
        })
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Failed to load example: {e}")


@router.post("/files/validate")
async def validate_session_file(data: LoadSessionFileInput):
    """
    Validate a session file without loading full data.

    Checks that all required fields are present and valid.
    """
    try:
        session = load_session_file(data.filepath)
        is_valid, errors, warnings = validate_session(session)

        return JSONResponse({
            "success": True,
            "filepath": data.filepath,
            "valid": is_valid,
            "errors": errors,
            "warnings": warnings,
            "summary": export_session_summary(session),
        })
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Validation failed: {e}")


@router.get("/files/summary/{example_name}")
async def get_example_summary(example_name: str):
    """
    Get a summary of an example session without full data.

    Useful for displaying in a list or preview.
    """
    examples_dir = Path(__file__).parent.parent.parent / "examples"
    filepath = examples_dir / f"{example_name}.sfppy.json"

    if not filepath.exists():
        raise HTTPException(status_code=404, detail=f"Example '{example_name}' not found")

    try:
        session = load_session_file(str(filepath))
        return JSONResponse({
            "success": True,
            "example_name": example_name,
            "summary": export_session_summary(session),
        })
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/files/convert-to-patankar")
async def convert_session_to_patankar(data: LoadSessionFileInput):
    """
    Convert a session file to Patankar solver inputs.

    Returns the data structures needed to run the simulation.
    """
    try:
        session = load_session_file(data.filepath)
        inputs = session_to_patankar_inputs(session)

        # Convert to serializable format
        result = {
            "geometry": {
                "type": str(type(inputs['geometry']).__name__),
                "volume_m3": float(inputs['volume']),
                "surface_m2": float(inputs['surface_area']),
            },
            "layers": [
                {
                    "index": l['index'],
                    "polymer": l['polymer'],
                    "thickness_m": l['thickness_m'],
                    "substance_id": l['substance_id'],
                    "C0": l['C0'],
                }
                for l in inputs['layers']
            ],
            "contact_steps": inputs['contact_steps'],
            "food": inputs['food'],
            "substances": {
                sid: {
                    "name": str(sobj) if sobj else None,
                    "available": sobj is not None,
                }
                for sid, sobj in inputs['substances'].items()
            },
        }

        return JSONResponse({
            "success": True,
            "filepath": data.filepath,
            "patankar_inputs": result,
        })
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except ImportError as e:
        raise HTTPException(status_code=500, detail=f"Patankar module error: {e}")
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Conversion failed: {e}")
