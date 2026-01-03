"""
patankar — SFPPy Core Migration Modeling Library
================================================

Safe Food Packaging in Python (SFPPy) - Mass transfer modeling
for food contact materials.

This module provides:
- Polymer layer definitions and transport properties
- Food simulant and contact condition classes
- Migration solvers (Patankar finite-volume method)
- PubChem integration for substance properties
- 3D geometry for packaging

@project: SFPPy — Safe Food Packaging in Python
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@agroparistech.fr
@license: MIT
@url: https://github.com/ovitrac/SFPPy

Example Usage
-------------
>>> from patankar import LDPE, ethanol, migrant, senspatankar
>>> # Or use submodule imports (always supported):
>>> from patankar.layer import LDPE, layer
>>> from patankar.food import ethanol50
>>> import patankar.layer as polymer  # module import
"""

__version__ = "1.50"
__author__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@agroparistech.fr"

# =============================================================================
# First, ensure submodules are loaded into sys.modules
# This allows both patterns to work:
#   - from patankar import LDPE          (direct import)
#   - import patankar.layer as polymer   (submodule import)
# =============================================================================
from . import layer as _layer_module
from . import food as _food_module
from . import migration as _migration_module
from . import geometry as _geometry_module
from . import loadpubchem as _loadpubchem_module
from . import useroverride as _useroverride_module

# =============================================================================
# Layer Module - Polymer definitions and transport properties
# =============================================================================
from patankar.layer import (
    # Core classes - note: 'layer' class imported with alias to avoid shadowing module
    layer as layer_class,
    layerLink,
    mesh,
    # Utility functions
    list_materials,
    material,
    resolve_material,
    help_layer,
    list_layer_subclasses,
    _toSI,
    check_units,
    format_scientific_latex,
    # Widgets (for Jupyter)
    create_polymer_dropdown,
    create_multi_layer_widget,
    # Polymer classes - Polyethylenes
    LDPE,
    HDPE,
    LLDPE,
    # Polymer classes - Polypropylenes
    PP,
    PPrubber,
    oPP,
    # Polymer classes - Styrenics
    PS,
    rPS,
    HIPS,
    rHIPS,
    SBS,
    # Polymer classes - PVC
    rigidPVC,
    plasticizedPVC,
    # Polymer classes - PET
    gPET,
    wPET,
    # Polymer classes - Others
    PMMA,
    PVAc,
)

# =============================================================================
# Food Module - Simulants and contact conditions
# =============================================================================
from patankar.food import (
    # Base classes
    foodphysics,
    foodlayer,
    foodproperty,
    # Contact type classes
    realfood,
    simulant,
    nofood,
    setoff,
    realcontact,
    testcontact,
    # Texture classes
    texture,
    solid,
    semisolid,
    liquid,
    perfectlymixed,
    # Chemical affinity classes
    chemicalaffinity,
    fat,
    aqueous,
    intermediate,
    # Temperature condition classes
    frozen,
    chilled,
    ambient,
    transportation,
    hotambient,
    # Simulant instances (EU regulation)
    water,
    water3aceticacid,
    ethanol,
    ethanol50,
    ethanol95,
    oliveoil,
    isooctane,
    tenax,
    # Utility functions
    help_food,
    list_food_classes,
    create_food_tree_widget,
)

# Also import specific contact scenarios if available
try:
    from patankar.food import hotfilled, stacked
except ImportError:
    pass  # Not all versions have these

# =============================================================================
# Migration Module - Solvers and result containers
# =============================================================================
from patankar.migration import (
    # Main solver
    senspatankar,
    # Result containers
    CFSimulationContainer,
    Cprofile,
    SensPatankarResult,
    SensPatankarResultCollection,
    # Printing/export functions
    print_figure,
    print_pdf,
    print_png,
    print_svg,
    # Utilities
    autoname,
    colormap,
    rgb,
    # Widgets (for Jupyter)
    create_simulation_widget,
    create_plotmigration_widget,
)

# =============================================================================
# Geometry Module - 3D packaging shapes
# =============================================================================
from patankar.geometry import (
    # Main class
    Packaging3D,
    # Base classes
    Shape3D,
    Connector,
    # Shape classes
    Cylinder,
    Cone,
    OpenCone,
    OpenCylinder1,
    OpenCylinder2,
    RectangularPrism,
    OpenPrism1,
    OpenPrism2,
    Sphere,
    Hemisphere,
    SquarePyramid,
    OpenSquare1,
    OpenSquare2,
    # Utility functions
    help_geometry,
    get_geometries_and_synonyms,
    get_all_shapes_info,
    create_packaging_widget,
)

# =============================================================================
# LoadPubChem Module - Substance database integration
# =============================================================================
from patankar.loadpubchem import (
    # Main classes
    migrant,
    migrantToxtree,
    CompoundIndex,
    # Utility functions
    get_default_index,
    polarity_index,
    is_java_available,
    get_java_version,
    create_substance_widget,
)

# =============================================================================
# UserOverride Module - Custom property overrides
# =============================================================================
from patankar.useroverride import useroverride

# =============================================================================
# Private utilities exposed for user convenience
# =============================================================================
from patankar.private.mstruct import (
    # Matlab-like structure classes
    struct,
    param,
    paramauto,
    # Path string class
    pstr,
)

# =============================================================================
# Convenience aliases
# =============================================================================
# Solver alias (common usage pattern)
solver = senspatankar

# Store alias (common usage pattern)
store = CFSimulationContainer

# Note: We intentionally do NOT expose 'layer' at the top level because
# it conflicts with the submodule name. Users should use:
#   - from patankar.layer import layer   (for the class)
#   - import patankar.layer as polymer   (for the module)
# The layer_class is available internally if needed.

# Alias for backward compatibility with survey module
# Note: _food_module is already loaded, expose it as dfood
dfood = _food_module

# =============================================================================
# Public API
# =============================================================================
__all__ = [
    # Version info
    "__version__",
    "__author__",
    "__email__",
    # Layer module
    # Note: 'layer' class not exported to avoid conflict with submodule
    # Use: from patankar.layer import layer
    "layerLink",
    "mesh",
    "list_materials",
    "material",
    "resolve_material",
    "help_layer",
    "_toSI",
    "check_units",
    "LDPE",
    "HDPE",
    "LLDPE",
    "PP",
    "PPrubber",
    "oPP",
    "PS",
    "rPS",
    "HIPS",
    "rHIPS",
    "SBS",
    "gPET",
    "wPET",
    "rigidPVC",
    "plasticizedPVC",
    "PMMA",
    "PVAc",
    # Food module
    "foodphysics",
    "foodlayer",
    "foodproperty",
    "realfood",
    "simulant",
    "nofood",
    "setoff",
    "texture",
    "solid",
    "semisolid",
    "liquid",
    "perfectlymixed",
    "fat",
    "aqueous",
    "intermediate",
    "frozen",
    "chilled",
    "ambient",
    "transportation",
    "hotambient",
    "water",
    "water3aceticacid",
    "ethanol",
    "ethanol50",
    "ethanol95",
    "oliveoil",
    "isooctane",
    "tenax",
    "help_food",
    # Note: 'food' not exported to avoid conflict with submodule
    # Use: import patankar.food as food
    "dfood",
    # Migration module
    "senspatankar",
    "solver",
    "CFSimulationContainer",
    "store",
    "Cprofile",
    "SensPatankarResult",
    "SensPatankarResultCollection",
    "print_figure",
    "print_pdf",
    "print_png",
    "print_svg",
    "autoname",
    "colormap",
    "rgb",
    # Geometry module
    "Packaging3D",
    "Shape3D",
    "Cylinder",
    "Sphere",
    "RectangularPrism",
    "help_geometry",
    # LoadPubChem module
    "migrant",
    "migrantToxtree",
    "CompoundIndex",
    "polarity_index",
    # UserOverride module
    "useroverride",
    # Private utilities (mstruct)
    "struct",
    "param",
    "paramauto",
    "pstr",
]
