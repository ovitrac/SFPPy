"""
Session File I/O for SFPPy Studio

Handles loading, saving, and validating .sfppy.json session files.
Includes conversion to/from Patankar solver objects.

@author: SFPPy Studio
@license: MIT
"""

import sys
from pathlib import Path
from typing import Optional, Dict, List, Any, Tuple
from datetime import datetime
import json

# Add parent paths
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from studio.app.models.session import (
    Session,
    Metadata,
    Geometry,
    GeometryShape,
    GeometryDimensions,
    Substance,
    SubstanceSource,
    SubstanceProperties,
    Layer,
    LayerSubstance,
    PolymerType,
    Food,
    FoodCategory,
    FoodTexture,
    FoodAffinity,
    ContactStep,
    ContactType,
    SimulationConfig,
    Results,
    SubstanceResult,
    ValueWithUnit,
    create_empty_session,
)
from studio.app.utils.units import convert_to_si, convert_length, convert_time, convert_temperature


# =============================================================================
# Session File Operations
# =============================================================================

def load_session_file(filepath: str) -> Session:
    """
    Load a session from a .sfppy.json file.

    Args:
        filepath: Path to the session file

    Returns:
        Session object

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file is not valid JSON or doesn't match schema
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Session file not found: {filepath}")

    return Session.load(str(path))


def save_session_file(session: Session, filepath: str, overwrite: bool = False) -> str:
    """
    Save a session to a .sfppy.json file.

    Args:
        session: Session object to save
        filepath: Target file path
        overwrite: Whether to overwrite existing file

    Returns:
        Path to saved file

    Raises:
        FileExistsError: If file exists and overwrite=False
    """
    path = Path(filepath)

    # Ensure .sfppy.json extension
    if not str(path).endswith('.sfppy.json'):
        path = path.with_suffix('.sfppy.json')

    if path.exists() and not overwrite:
        raise FileExistsError(f"File already exists: {path}")

    # Update modified timestamp
    session.metadata.modified = datetime.utcnow()

    session.save(str(path))
    return str(path)


def list_session_files(directory: str, pattern: str = "*.sfppy.json") -> List[Dict[str, Any]]:
    """
    List all session files in a directory.

    Args:
        directory: Directory to search
        pattern: Glob pattern for matching files

    Returns:
        List of dicts with file info (path, name, modified, size)
    """
    dir_path = Path(directory)
    if not dir_path.exists():
        return []

    files = []
    for path in dir_path.glob(pattern):
        try:
            session = Session.load(str(path))
            files.append({
                "path": str(path),
                "filename": path.name,
                "name": session.metadata.name,
                "description": session.metadata.description,
                "created": session.metadata.created.isoformat() if session.metadata.created else None,
                "modified": session.metadata.modified.isoformat() if session.metadata.modified else None,
                "substances_count": len(session.substances),
                "layers_count": len(session.layers),
                "steps_count": len(session.contact_steps),
                "size_bytes": path.stat().st_size,
            })
        except Exception as e:
            files.append({
                "path": str(path),
                "filename": path.name,
                "error": str(e),
            })

    return files


# =============================================================================
# Session Validation
# =============================================================================

def validate_session(session: Session) -> Tuple[bool, List[str], List[str]]:
    """
    Validate a session is ready for simulation.

    Args:
        session: Session to validate

    Returns:
        Tuple of (is_valid, errors, warnings)
    """
    errors = []
    warnings = []

    # Check layers
    if not session.layers:
        errors.append("At least one layer is required")
    else:
        for layer in session.layers:
            if layer.thickness.value <= 0:
                errors.append(f"Layer {layer.index}: thickness must be positive")

    # Check contact steps
    if not session.contact_steps:
        errors.append("At least one contact step is required")
    else:
        for step in session.contact_steps:
            if step.duration.value <= 0:
                errors.append(f"Step {step.index}: duration must be positive")

    # Check substances are assigned
    has_any_assignment = False
    for layer in session.layers:
        if layer.substances:
            for sub in layer.substances:
                if sub.C0.value > 0:
                    has_any_assignment = True
                    break

    if not has_any_assignment:
        warnings.append("No substance with C0 > 0 assigned to any layer")

    # Check geometry
    if session.geometry.shape == GeometryShape.CYLINDER:
        if not session.geometry.dimensions.radius:
            errors.append("Cylinder geometry requires radius")
        if not session.geometry.dimensions.height and not session.geometry.dimensions.length:
            errors.append("Cylinder geometry requires height or length")
    elif session.geometry.shape == GeometryShape.BOX_CONTAINER:
        if not session.geometry.dimensions.length:
            errors.append("Box geometry requires length")
        if not session.geometry.dimensions.width:
            errors.append("Box geometry requires width")
        if not session.geometry.dimensions.height:
            errors.append("Box geometry requires height")

    is_valid = len(errors) == 0
    return is_valid, errors, warnings


# =============================================================================
# Conversion to Patankar Objects
# =============================================================================

def session_to_patankar_inputs(session: Session) -> Dict[str, Any]:
    """
    Convert a Session to Patankar solver input objects.

    Returns a dictionary with keys:
    - 'geometry': Packaging3D object
    - 'layers': List of layer objects
    - 'food': Food object
    - 'substances': Dict of migrant objects
    - 'contact_steps': List of contact condition dicts

    Args:
        session: Session object

    Returns:
        Dictionary of Patankar-compatible objects
    """
    try:
        from patankar.geometry import Packaging3D
        from patankar.loadpubchem import migrant
        import patankar.layer as polymer_module
        import patankar.food as food_module
    except ImportError as e:
        raise ImportError(f"Patankar module not available: {e}")

    result = {}

    # --- Geometry ---
    geom = session.geometry
    if geom.shape == GeometryShape.CYLINDER:
        radius = geom.dimensions.radius.value if geom.dimensions.radius else 30
        radius_unit = geom.dimensions.radius.unit if geom.dimensions.radius else "mm"
        height = geom.dimensions.height.value if geom.dimensions.height else (
            geom.dimensions.length.value if geom.dimensions.length else 100
        )
        height_unit = geom.dimensions.height.unit if geom.dimensions.height else (
            geom.dimensions.length.unit if geom.dimensions.length else "mm"
        )
        result['geometry'] = Packaging3D(
            'Cylinder',
            radius=(radius, radius_unit),
            height=(height, height_unit)
        )
    elif geom.shape == GeometryShape.BOX_CONTAINER:
        result['geometry'] = Packaging3D(
            'box_container',
            length=(geom.dimensions.length.value, geom.dimensions.length.unit),
            width=(geom.dimensions.width.value, geom.dimensions.width.unit),
            height=(geom.dimensions.height.value, geom.dimensions.height.unit)
        )
    else:
        # Default to cylinder
        result['geometry'] = Packaging3D('Cylinder', radius=(30, "mm"), height=(100, "mm"))

    # Get volume and surface area
    result['volume'], result['surface_area'] = result['geometry'].get_volume_and_area()

    # --- Substances ---
    substances = {}
    for sub in session.substances:
        if sub.source == SubstanceSource.PUBCHEM and sub.lookup_name:
            try:
                m = migrant(sub.lookup_name)
                substances[sub.id] = m
            except Exception:
                # Create placeholder if lookup fails
                substances[sub.id] = None
        else:
            substances[sub.id] = None
    result['substances'] = substances

    # --- Layers ---
    # Map polymer types to Patankar classes
    polymer_map = {
        PolymerType.LDPE: polymer_module.LDPE,
        PolymerType.LLDPE: polymer_module.LLDPE,
        PolymerType.HDPE: polymer_module.HDPE,
        PolymerType.PP: polymer_module.PP,
        PolymerType.gPET: polymer_module.gPET,
        PolymerType.wPET: polymer_module.wPET,
        PolymerType.rPET: polymer_module.rPET,
        PolymerType.PS: polymer_module.PS,
        PolymerType.PMMA: polymer_module.PMMA,
        PolymerType.PA6: polymer_module.PA6,
        PolymerType.PA66: polymer_module.PA66,
        PolymerType.PBT: polymer_module.PBT,
        PolymerType.PEN: polymer_module.PEN,
        PolymerType.HIPS: polymer_module.HIPS,
        PolymerType.SBS: polymer_module.SBS,
    }

    layers = []
    for layer in sorted(session.layers, key=lambda x: x.index):
        polymer_class = polymer_map.get(layer.polymer, polymer_module.LDPE)
        thickness_m = convert_length(layer.thickness.value, layer.thickness.unit, "m")

        # Get substance assignment for this layer
        sub_assignment = layer.substances[0] if layer.substances else None
        sub_obj = None
        c0 = 0

        if sub_assignment:
            sub_obj = substances.get(sub_assignment.substance_id)
            c0 = sub_assignment.C0.value

        # Create layer
        layer_obj = polymer_class(
            l=(thickness_m, "m"),
            migrant=sub_obj,
            C0=c0
        )
        layers.append({
            'index': layer.index,
            'object': layer_obj,
            'polymer': layer.polymer.value,
            'thickness_m': thickness_m,
            'substance_id': sub_assignment.substance_id if sub_assignment else None,
            'C0': c0,
        })

    result['layers'] = layers

    # --- Contact Steps ---
    contact_steps = []
    for step in sorted(session.contact_steps, key=lambda x: x.index):
        temp_C = convert_temperature(step.temperature.value, step.temperature.unit, "degC")
        duration_s = convert_time(step.duration.value, step.duration.unit, "s")

        contact_steps.append({
            'index': step.index,
            'name': step.name,
            'type': step.type.value,
            'temperature_C': temp_C,
            'duration_s': duration_s,
            'with_food_contact': step.with_food_contact,
        })

    result['contact_steps'] = contact_steps

    # --- Food ---
    result['food'] = {
        'category': session.food.category.value,
        'texture': session.food.texture.value,
        'affinity': session.food.affinity.value,
        'simulant': session.food.simulant,
        'name': session.food.name,
    }

    return result


# =============================================================================
# Results Storage
# =============================================================================

def store_results_in_session(
    session: Session,
    results: Dict[str, Any],
    elapsed_seconds: float
) -> Session:
    """
    Store simulation results in the session.

    Args:
        session: Session object
        results: Dictionary with simulation results
        elapsed_seconds: Computation time

    Returns:
        Updated session with results
    """
    substance_results = []

    for sub_id, sub_result in results.get('substances', {}).items():
        # Find SML for this substance
        sub_def = next((s for s in session.substances if s.id == sub_id), None)
        sml = sub_def.SML if sub_def else None

        cf_target = sub_result.get('CF_target', 0)
        cf_equilibrium = sub_result.get('CF_equilibrium', cf_target)

        compliant = None
        margin = None
        if sml and sml > 0:
            compliant = cf_target <= sml
            margin = 100 * (sml - cf_target) / sml

        substance_results.append(SubstanceResult(
            substance_id=sub_id,
            CF_at_tcontact=cf_target,
            CF_equilibrium=cf_equilibrium,
            SML=sml,
            compliant=compliant,
            margin_percent=margin,
        ))

    session.results = Results(
        computed_at=datetime.utcnow(),
        elapsed_seconds=elapsed_seconds,
        substances=substance_results,
    )

    return session


# =============================================================================
# Export Functions
# =============================================================================

def export_session_summary(session: Session) -> Dict[str, Any]:
    """
    Export a session summary for display/reporting.

    Args:
        session: Session object

    Returns:
        Dictionary with summary information
    """
    # Calculate total contact time
    total_duration_s = sum(
        convert_time(step.duration.value, step.duration.unit, "s")
        for step in session.contact_steps
    )
    total_duration_days = total_duration_s / 86400

    # Calculate total thickness
    total_thickness_um = sum(
        convert_length(layer.thickness.value, layer.thickness.unit, "um")
        for layer in session.layers
    )

    # Get substance names
    substance_names = [
        s.properties.name if s.properties else s.lookup_name or s.id
        for s in session.substances
    ]

    return {
        "name": session.metadata.name,
        "description": session.metadata.description,
        "geometry": {
            "shape": session.geometry.shape.value,
        },
        "layers": {
            "count": len(session.layers),
            "total_thickness_um": total_thickness_um,
            "polymers": [l.polymer.value for l in session.layers],
        },
        "substances": {
            "count": len(session.substances),
            "names": substance_names,
        },
        "contact_steps": {
            "count": len(session.contact_steps),
            "total_duration_days": total_duration_days,
            "types": [s.type.value for s in session.contact_steps],
        },
        "food": {
            "category": session.food.category.value,
            "texture": session.food.texture.value,
            "affinity": session.food.affinity.value,
        },
        "has_results": session.results is not None,
    }
