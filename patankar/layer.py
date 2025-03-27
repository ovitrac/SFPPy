#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
SFPPy Module: Layer (Packaging Materials)
===============================================================================
Defines **packaging materials** as 1D layers. Supports:
- **Multilayer assembly (`layer1 + layer2`)**
- **Mass transfer modeling (`layer >> food`)**
- **Automatic meshing for finite-volume solvers**

**Main Components:**
- **Base Class: `layer`** (Defines all packaging materials)
    - Properties: `D` (diffusivity), `k` (partition coefficient), `l` (thickness)
    - Supports **+ (stacking)** and **splitting** layers
    - Propagates contact temperature from `food.py`
- **Predefined Materials (Subclasses)**:
    - `LDPE`, `PP`, `PET`, `Cardboard`, `Ink`
- **Dynamic Property Models:**
    - `Dmodel()`, `kmodel()`: Call `property.py` to predict diffusion and partitioning

**Integration with SFPPy Modules:**
- Used in `migration.py` to define the **left-side boundary**.
- Retrieves chemical properties from `loadpubchem.py`.
- Works with `food.py` to model **food-contact** interactions.

Example:
```python
from patankar.layer import LDPE
A = LDPE(l=50e-6, D=1e-14)
```


===============================================================================
Details
===============================================================================
Layer builder for patankar package

All materials are represented as layers and be combined, merged with mathematical
operations such as +. The general object general object is of class layer.

Specific materials with known properties have been derived: LDPE(),HDPE(),PP()...air()

List of implemented materials:

    | Class Name              | Type     | Material                        | Code    |
    |-------------------------|----------|---------------------------------|---------|
    | AdhesiveAcrylate        | adhesive | acrylate adhesive               | Acryl   |
    | AdhesiveEVA             | adhesive | EVA adhesive                    | EVA     |
    | AdhesiveNaturalRubber   | adhesive | natural rubber adhesive         | rubber  |
    | AdhesivePU              | adhesive | polyurethane adhesive           | PU      |
    | AdhesivePVAC            | adhesive | PVAc adhesive                   | PVAc    |
    | AdhesiveSyntheticRubber | adhesive | synthetic rubber adhesive       | sRubber |
    | AdhesiveVAE             | adhesive | VAE adhesive                    | VAE     |
    | Cardboard               | paper    | cardboard                       | board   |
    | HDPE                    | polymer  | high-density polyethylene       | HDPE    |
    | HIPS                    | polymer  | high-impact polystyrene         | HIPS    |
    | LDPE                    | polymer  | low-density polyethylene        | LDPE    |
    | LLDPE                   | polymer  | linear low-density polyethylene | LLDPE   |
    | PA6                     | polymer  | polyamide 6                     | PA6     |
    | PA66                    | polymer  | polyamide 6,6                   | PA6,6   |
    | SBS                     | polymer  | styrene-based polymer SBS       | SBS     |
    | PBT                     | polymer  | polybutylene terephthalate      | PBT     |
    | PEN                     | polymer  | polyethylene naphthalate        | PEN     |
    | PP                      | polymer  | isotactic polypropylene         | PP      |
    | PPrubber                | polymer  | atactic polypropylene           | aPP     |
    | PS                      | polymer  | polystyrene                     | PS      |
    | Paper                   | paper    | paper                           | paper   |
    | air                     | air      | ideal gas                       | gas     |
    | gPET                    | polymer  | glassy PET                      | PET     |
    | oPP                     | polymer  | bioriented polypropylene        | oPP     |
    | plasticizedPVC          | polymer  | plasticized PVC                 | pPVC    |
    | rPET                    | polymer  | rubbery PET                     | rPET    |
    | rigidPVC                | polymer  | rigid PVC                       | PVC     |


Mass transfer within each layer are governed by a diffusion coefficient, a Henri-like coefficient
enabling to describe the partitioning between layers. All materials are automatically meshed using
a modified finite volume technique exact at steady state and offering good accuracy in non-steady
conditions.

A temperature and substance can be assigned to layers.


@version: 1.40
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2022-02-21
@rev. 2025-03-26

"""

# ---- History ----
# Created on Tue Jan 18 09:14:34 2022
# 2025-01-19 RC
# 2022-01-20 full indexing and simplification
# 2022-01-21 add split()
# 2022-01-22 add child classes for common polymers
# 2022-01-23 full implementation of units
# 2022-01-26 mesh() method generating mesh objects
# 2022-02-21 add compatibility with migration


oneline = "Build multilayer objects"

docstr = """
Build layer(s) for SENSPATANKAR

Example of caller:
    from patankar.layer import layer
    A=layer(D=1e-14,l=50e-6)
    A

"""


# Package Dependencies
# ====================
# <--  generic packages  -->
import sys
import inspect
import textwrap
import numpy as np
from copy import deepcopy as duplicate
# <--  local packages  -->
if 'SIbase' not in dir(): # avoid loading it twice
    from patankar.private.pint import UnitRegistry as SIbase
    from patankar.private.pint import set_application_registry as fixSIbase
if 'migrant' not in dir():
    from patankar.loadpubchem import migrant
from patankar.useroverride import useroverride # useroverride is already an instance (not a class)


__all__ = ['AdhesiveAcrylate', 'AdhesiveEVA', 'AdhesiveNaturalRubber', 'AdhesivePU', 'AdhesivePVAC', 'AdhesiveSyntheticRubber', 'AdhesiveVAE', 'Cardboard', 'HDPE', 'HIPS', 'LDPE', 'LLDPE', 'PA6', 'PA66', 'PBT', 'PEN', 'PMMA', 'PP', 'PPrubber', 'PS', 'PVAc', 'Paper', 'R', 'RT0K', 'SBS', 'SI', 'SIbase', 'T0K', 'air', 'check_units', 'create_multi_layer_widget', 'create_polymer_dropdown', 'fixSIbase', 'format_scientific_latex', 'gPET', 'help_layer', 'iRT0K', 'layer', 'layerLink', 'list_layer_subclasses', 'list_materials', 'material', 'mesh', 'migrant', 'oPP', 'plasticizedPVC', 'qSI', 'rHIPS', 'rPET', 'rPS', 'resolve_material', 'rigidPVC', 'toSI', 'useroverride', 'wPET']

__project__ = "SFPPy"
__author__ = "Olivier Vitrac"
__copyright__ = "Copyright 2022"
__credits__ = ["Olivier Vitrac"]
__license__ = "MIT"
__maintainer__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@agroparistech.fr"
__version__ = "1.40"

# %% Private functions and classes

# Unified registry of materials
# Each key is the canonical class name (as used in layer.py) and its value is a dictionary with:
#   - 'description': the full display name
#   - 'type': one of "polymer", "adhesive", or "other"
#   - 'synonyms': a list of alternative names that can be used for import purposes

mainregistry = {
    # Polymers
    'HDPE': {
        'description': 'High-Density Polyethylene',
        'type': 'polymer',
        'synonyms': ["PEHD","high-density polyethylene"]
    },
    'HIPS': {
        'description': 'High-Impact Polystyrene',
        'type': 'polymer',
        'synonyms': ["high-impact polystyrene"]
    },
    'LDPE': {
        'description': 'Low-Density Polyethylene',
        'type': 'polymer',
        'synonyms': ["PEBD","PE","low-desnity polyethylne"]
    },
    'LLDPE': {
        'description': 'Linear Low-Density Polyethylene',
        'type': 'polymer',
        'synonyms': ["PEBLD","PELBD", "linear low-Density polyethylene"]
    },
    'PA6': {
        'description': 'Polyamide 6',
        'type': 'polymer',
        'synonyms': ["PA 6","Polyamide 6"]
    },
    'PA66': {
        'description': 'Polyamide 6,6',
        'type': 'polymer',
        'synonyms': ["PA 6,6","polyamide 6,6"]
    },
    'PBT': {
        'description': 'Polybutylene Terephthalate',
        'type': 'polymer',
        'synonyms': []
    },
    'PEN': {
        'description': 'Polyethylene Naphthalate',
        'type': 'polymer',
        'synonyms': ["polyethylene naphthalate"]
    },
    'PMMA': {
        'description': 'Polymethyl Methacrylate',
        'type': 'polymer',
        'synonyms': ["polymethyl methacrylate"]
    },
    'PP': {
        'description': 'Polypropylene',
        'type': 'polymer',
        'synonyms': ["polypropylene"]
    },
    'PPrubber': {
        'description': 'Atactic Polypropylene (Rubbery)',
        'type': 'polymer',
        'synonyms': ["aPP","atactic PP","atactic polypropylene"]
    },
    'PS': {
        'description': 'Polystyrene',
        'type': 'polymer',
        'synonyms': ["GPPS","polystyrene","general purpose polystyrene"]
    },
    'PVAc': {
        'description': 'Polyvinyl Acetate',
        'type': 'polymer',
        'synonyms': ["PVAC"]
    },
    'SBS': {
        'description': 'Styrene-Butadiene-Styrene',
        'type': 'polymer',
        'synonyms': ["styrene-butadiene-styrene"]
    },
    'gPET': {
        'description': 'Glassy Polyethylene Terephthalate',
        'type': 'polymer',
        'synonyms': ['PET', 'PETE', 'polyethylene terephthalate']
    },
    'oPP': {
        'description': 'Oriented Polypropylene',
        'type': 'polymer',
        'synonyms': ["oriented polypropylene"]
    },
    'plasticizedPVC': {
        'description': 'Plasticized Polyvinyl Chloride',
        'type': 'polymer',
        'synonyms': ["plasticized polyvinyl chloride","PVC cling film"]
    },
    'rHIPS': {
        'description': 'Rubbery High-Impact Polystyrene',
        'type': 'polymer',
        'synonyms': ["rubbery HIPS"]
    },
    'rPET': {
        'description': 'Rubbery Polyethylene Terephthalate',
        'type': 'polymer',
        'synonyms': ['PET above Tg', 'PETE above Tg', 'PET T>Tg', 'PETE T>Tg']
    },
    'rPS': {
        'description': 'Rubbery Polystyrene',
        'type': 'polymer',
        'synonyms': []
    },
    'rigidPVC': {
        'description': 'Rigid Polyvinyl Chloride',
        'type': 'polymer',
        'synonyms': ["rigid PVC","PVC","rigid polyvinyl chloride"]
    },
    'wPET': {
        'description': 'Wet (Plasticized) Polyethylene Terephthalate',
        'type': 'polymer',
        'synonyms': ['plasticized PET', 'plasticized PETE','wet PET']
    },

    # Adhesives
    'AdhesiveAcrylate': {
        'description': 'Acrylate Adhesive',
        'type': 'adhesive',
        'synonyms': ["PMMA adhesive"]
    },
    'AdhesiveEVA': {
        'description': 'EVA Adhesive',
        'type': 'adhesive',
        'synonyms': ["EVA adhesive"]
    },
    'AdhesiveNaturalRubber': {
        'description': 'Natural Rubber Adhesive',
        'type': 'adhesive',
        'synonyms': ["natural rubber adhesive"]
    },
    'AdhesivePU': {
        'description': 'Polyurethane Adhesive',
        'type': 'adhesive',
        'synonyms': ["PU adhesive","polyurethane adhesive"]
    },
    'AdhesivePVAC': {
        'description': 'PVAc Adhesive',
        'type': 'adhesive',
        'synonyms': ["adhesive PVAc","adhesive PVAC"]
    },
    'AdhesiveSyntheticRubber': {
        'description': 'Synthetic Rubber Adhesive',
        'type': 'adhesive',
        'synonyms': ["synthetic rubber adhesive"]
    },
    'AdhesiveVAE': {
        'description': 'VAE Adhesive',
        'type': 'adhesive',
        'synonyms': ["VAE adhesive"]
    },

    # Other materials
    'Cardboard': {
        'description': 'Cardboard',
        'type': 'other',
        'synonyms': ["board"]
    },
    'Paper': {
        'description': 'Paper',
        'type': 'other',
        'synonyms': ["paper"]
    },
    'air': {
        'description': 'Air',
        'type': 'other',
        'synonyms': ["gaz"]
    }
}

# Listing all available materials with details
def list_materials():
    """
    Lists all available materials with their descriptions, types, and synonyms
    as a Markdown table with adjusted column widths.
    """
    headers = ["Material Key", "Description", "Type", "Synonyms"]
    rows = []
    for key, info in mainregistry.items():
        synonyms = ", ".join(info["synonyms"]) if info["synonyms"] else "None"
        rows.append([key, info["description"], info["type"], synonyms])
    # Compute the maximum width for each column based on headers and rows
    col_widths = [
        max(len(headers[i]), *(len(row[i]) for row in rows))
        for i in range(len(headers))
    ]
    # Helper to format a row given the computed widths
    def format_row(row):
        return "| " + " | ".join(row[i].ljust(col_widths[i]) for i in range(len(row))) + " |"
    # Build the Markdown table rows
    table = [format_row(headers)]
    # Separator row: use dashes matching the column widths
    separator = "| " + " | ".join("-" * col_widths[i] for i in range(len(headers))) + " |"
    table.append(separator)
    # Data rows
    for row in rows:
        table.append(format_row(row))

    # Print the complete Markdown table
    print("\n".join(table))


# Import mechanism based on synonyms
def resolve_material(name,returnclass=True):
    """
    Resolves a material name or any of its synonyms to the canonical registry key.

    Parameters:
      name (str): The material name or synonym.

    Returns:
      key (str): The canonical key from the registry.

    Raises:
      KeyError: If no matching material is found.
    """
    name_lower = name.lower()
    for key, info in mainregistry.items():
        # Check if the provided name matches the canonical key (case-insensitive)
        if key.lower() == name_lower:
            return key
        # Check if the provided name matches any synonym (case-insensitive)
        for syn in info['synonyms']:
            if syn.lower() == name_lower:
                return key
    print("List of available materials.\nIf you not find yours, call layer directly for a customized one.")
    list_materials()
    raise KeyError(f"Material '{name}' not found in the registry.")


def material(name):
    """
    Import surrogate that returns the material class corresponding to the given name.

    Usage:
        from patankar import material
        mymaterial = material("anyname")  # 'anyname' is resolved to a canonical key,
                                          # and the corresponding class (or the base 'layer' if applicable)
                                          # from this module is returned.

    If the provided name is "layer" (case-insensitive), the base layer class is returned.

    Parameters:
      name (str): The material name or any synonym.

    Returns:
      A material class that can be instantiated (e.g., mymaterial(l=...)).

    Raises:
      KeyError: If the provided name (or its synonym) does not match any material in the registry.
      ImportError: If the material class cannot be found in the module.

    Example:
        A=material("PEBD")(l=(10,"µm")) # instantiate LDPE with l=100 µm
    """
    # If the user requests the base class "layer", return it directly.
    if name.lower() == "layer":
        return layer
    # Resolve the canonical key using the existing resolve_material function.
    canonical_key = resolve_material(name)
    try:
        # Since material is defined in this module, use globals() to fetch the class.
        mat_class = globals()[canonical_key]
    except KeyError:
        raise ImportError(f"Material class '{canonical_key}' not found in the module.")
    return mat_class



# Initialize unit conversion (intensive initialization with old Python versions)
# NB: degC and kelvin must be used for temperature
# conversion as obj,valueSI,unitSI = toSI(qSI(numvalue,"unit"))
# conversion as obj,valueSI,unitSI = toSI(qSI("value unit"))
def toSI(q): q=q.to_base_units(); return q,q.m,str(q.u)
NoUnits = 'a.u.'     # value for arbitrary unit
UnknownUnits = 'N/A' # no non indentified units
if ("SI" not in locals()) or ("qSI" not in locals()):
    SI = SIbase()      # unit engine
    fixSIbase(SI)      # keep the same instance between calls
    qSI = SI.Quantity  # main unit consersion method from string
    # constants (usable in layer object methods)
    # define R,T0K,R*T0K,1/(R*T0K) with there SI units
    constants = {}
    R,constants["R"],constants["Runit"] = toSI(qSI(1,'avogadro_number*boltzmann_constant'))
    T0K,constants["T0K"],constants["T0Kunit"] = toSI(qSI(0,'degC'))
    RT0K,constants["RT0K"],constants["RT0Kunit"] = toSI(R*T0K)
    iRT0K,constants["iRT0K"],constants["iRT0Kunit"] = toSI(1/RT0K)

# Concise data validator with unit convertor to SI
# To prevent many issues with temperature and to adhere to 2024 golden standard in layer
# defaulttempUnits has been set back to "degC" from "K".
def check_units(value,ProvidedUnits=None,ExpectedUnits=None,defaulttempUnits="degC"):
    """ check numeric inputs and convert them to SI units """
    # by convention, NumPy arrays and None are return unchanged (prevent nesting)
    if isinstance(value,np.ndarray) or value is None:
        return value,UnknownUnits
    if isinstance(value,tuple):
        if len(value) != 2:
            raise ValueError('value should be a tuple: (value,"unit"')
        ProvidedUnits = value[1]
        value = value[0]
    if isinstance(value,list): # the function is vectorized
        value = np.array(value)
    if {"degC", "K"} & {ProvidedUnits, ExpectedUnits}: # the value is a temperature
        ExpectedUnits = defaulttempUnits if ExpectedUnits is None else ExpectedUnits
        ProvidedUnits = ExpectedUnits if ProvidedUnits is None else ProvidedUnits
        if ProvidedUnits=="degC" and ExpectedUnits=="K":
            value += constants["T0K"]
        elif ProvidedUnits=="K" and ExpectedUnits=="degC":
            value -= constants["T0K"]
        return np.array([value]),ExpectedUnits
    else: # the value is not a temperature
        ExpectedUnits = NoUnits if ExpectedUnits is None else ExpectedUnits
        if (ProvidedUnits==ExpectedUnits) or (ProvidedUnits==NoUnits) or (ExpectedUnits==None):
            conversion =1               # no conversion needed
            units = ExpectedUnits if ExpectedUnits is not None else NoUnits
        else:
            q0,conversion,units = toSI(qSI(1,ProvidedUnits))
        return np.array([value*conversion]),units

# _toSI: function helper for the enduser outside layer
def _toSI(value=None):
    '''return an SI value from (value,"unit")'''
    if not isinstance(value,tuple) or len(value)!=2 \
        or not isinstance(value[0],(float,int,list,np.ndarray)) \
            or  not isinstance(value[1],str):
        raise ValueError('value must be (currentvalue,"unit") - for example: (10,"days")')
    return check_units(value)[0]


# formatsci equivalent
def format_scientific_latex(value, numdigits=4, units=None, prefix="",mathmode="$"):
    """
    Formats a number in scientific notation only when necessary, using LaTeX.

    Parameters:
    -----------
    value : float
        The number to format.
    numdigits : int, optional (default=4)
        Number of significant digits for formatting.
    units : str, optional (default=None)
        LaTeX representation of units. If None, no units are added.
    prefix: str, optional (default="")
    mathmode: str, optional (default="$")

    Returns:
    --------
    str
        The formatted number in standard or LaTeX scientific notation.

    Examples:
    ---------
    >>> format_scientific_latex(1e-12)
    '$10^{-12}$'

    >>> format_scientific_latex(1.5e-3)
    '0.0015'

    >>> format_scientific_latex(1.3e10)
    '$1.3 \\cdot 10^{10}$'

    >>> format_scientific_latex(0.00341)
    '0.00341'

    >>> format_scientific_latex(3.41e-6)
    '$3.41 \\cdot 10^{-6}$'
    """

    if value == 0:
        return "$0$" if units is None else rf"$0 \, {units}$"
    # Get formatted number using Matlab-like %g behavior
    formatted = f"{value:.{numdigits}g}"
    # If the formatting results in an `e` notation, convert to LaTeX
    if "e" in formatted or "E" in formatted:
        coefficient, exponent = formatted.split("e")
        exponent = int(exponent)  # Convert exponent to integer
        # Remove trailing zeros in coefficient
        coefficient = coefficient.rstrip("0").rstrip(".")  # Ensures "1.00" -> "1"
        # LaTeX scientific format
        sci_notation = rf"{prefix}{coefficient} \cdot 10^{{{exponent}}}"
        return sci_notation if units is None else rf"{mathmode}{sci_notation} \, {units}{mathmode}"
    # Otherwise, return standard notation
    return formatted if units is None else rf"{mathmode}{prefix}{formatted} \, {units}{mathmode}"




# helper function to list all classes
def list_layer_subclasses():
    """
    Lists all classes in this module that derive from 'layer',
    along with their layertype and layermaterial properties.

    Returns:
        list of tuples (classname, layertype, layermaterial)
    """
    subclasses_info = []
    current_module = sys.modules[__name__]  # This refers to layer.py itself
    for name, obj in inspect.getmembers(current_module, inspect.isclass):
        # Make sure 'obj' is actually a subclass of layer (and not 'layer' itself)
        if obj is not layer and issubclass(obj, layer):
            try:
                # Instantiate with default parameters so that .layertype / .layermaterial are accessible
                instance = obj()
                subclasses_info.append(
                    {"classname":name,
                     "type":instance._type[0],
                     "material":instance._material[0],
                     "code":instance._code[0]}
                )
            except TypeError as e:
                # Log error and rethrow for debugging
                print(f"⚠️ Error: Could not instantiate class '{name}'. Check its constructor.")
                print(f"🔍 Exception: {e}")
                raise  # Rethrow the error with full traceback
    return subclasses_info


# general help for layer
def help_layer():
    """
    Print all subclasses with their type/material info in a Markdown table with dynamic column widths.
    """
    derived = list_layer_subclasses()
    # Extract table content
    headers = ["Class Name", "Type", "Material", "Code"]
    rows = [[item["classname"], item["type"], item["material"], item["code"]] for item in derived]
    # Compute column widths based on content
    col_widths = [max(len(str(cell)) for cell in col) for col in zip(headers, *rows)]
    # Formatting row template
    row_format = "| " + " | ".join(f"{{:<{w}}}" for w in col_widths) + " |"
    # Print header
    print(row_format.format(*headers))
    print("|-" + "-|-".join("-" * w for w in col_widths) + "-|")

    # Print table rows
    for row in rows:
        print(row_format.format(*row))


# %% Widgets

#ipywidgets: polymer selection
def create_polymer_dropdown(default_value=None):
    """
    Creates and returns a dropdown widget for selecting a polymer.
    The dropdown options are built from the mainregistry entries with type 'polymer',
    formatted as "Key: Description". If ipywidgets is not installed, a dummy widget
    is returned that mimics the basic interface.

    Parameters:
      default_value (str, optional): The default polymer key to select.
                                     If not provided, the first polymer key from the registry is used.

    Returns:
      A widget (ipywidgets.Dropdown or a fallback DummyDropdown) with the polymer options.
    """
    # Build polymer options as a list of (display, value) tuples
    polymer_options = []
    for key, info in mainregistry.items():
        if info["type"] == "polymer":
            display_text = f"{key}: {info['description']}"
            polymer_options.append((display_text, key))

    if not polymer_options:
        raise ValueError("No polymers found in the registry.")

    # Set default value if not provided or if not in options
    valid_keys = [opt[1] for opt in polymer_options]
    if default_value is None or default_value not in valid_keys:
        default_value = valid_keys[0]

    # Try to import ipywidgets; if unavailable, define a fallback DummyDropdown.
    try:
        import ipywidgets as widgets
    except ImportError:
        widgets = None

    if widgets:
        dropdown = widgets.Dropdown(
            options=polymer_options,
            value=default_value,
            description='Polymer:',
            layout=widgets.Layout(width='50%')
        )
    else:
        # Fallback dummy widget mimicking ipywidgets.Dropdown
        class DummyDropdown:
            def __init__(self, options, value, description, layout):
                self.options = options
                self.value = value
                self.description = description
                self.layout = layout

            def __repr__(self):
                return (f"DummyDropdown(value={self.value}, "
                        f"options={[opt[1] for opt in self.options]}, "
                        f"description='{self.description}', layout={self.layout})")
        dropdown = DummyDropdown(
            options=polymer_options,
            value=default_value,
            description='Polymer:',
            layout={'width': '50%'}
        )

    return dropdown



def create_multi_layer_widget(default_polymer="LDPE", default_thickness_value=100, default_thickness_unit="um", default_c0=100):
    """
    Creates a widget interface to define multiple layers (1 to 10) with the following parameters for each layer:
      - Layer Name (unique identifier)
      - Polymer (dropdown with "Key: Description")
      - Thickness (a numeric value and a unit)
      - Initial Concentration C0 (in arbitrary units, a.u.)

    The interface includes:
      - An IntSlider to choose the number of layers.
      - Navigation arrows (Previous/Next) with a label showing "Layer X of Y".
      - Input fields for each layer's parameters.
      - A button to instantiate all layers. When clicked, each layer is instantiated as:
            material(layer_def["polymer"])(l=(thickness, unit), D=1e-14, C0=(C0, "a.u."))
        and stored in a global dictionary (builtins.mylayers) keyed by the layer name.

    The user may add or remove layers by changing the number and editing the layer names.

    Returns:
      An ipywidgets.VBox instance containing the full UI.

    Raises:
      ImportError: If ipywidgets/IPython are not available.
    """
    try:
        import ipywidgets as widgets
        from IPython.display import display
    except ImportError as e:
        raise ImportError("ipywidgets and IPython are required for the widget interface.") from e

    import builtins
    # Ensure a global dictionary for storing created layers exists
    if not hasattr(builtins, "mylayers"):
        builtins.mylayers = {}
    # Instead of mymaterial, we now use mymaterials:
    if not hasattr(builtins, "mymaterials"):
        builtins.mymaterials = {}
    global mylayers, mymaterials
    mylayers = builtins.mylayers
    mymaterials = builtins.mymaterials
    # flag for preheated gui interface (widgets should be initialized manually, instead of being empty)
    _preheatedGUI_ = hasattr(builtins, "_PREHEATED_") and getattr(builtins, "_PREHEATED_") is True

    # Widget to select number of layers
    num_layers_slider = widgets.IntSlider(
        value=1,
        min=1,
        max=10,
        step=1,
        description='Num Layers:',
        continuous_update=False,
        layout=widgets.Layout(width='50%')
    )

    # Navigation label (e.g., "Layer 1 of 3")
    nav_label = widgets.Label(value="Layer 1 of 1")

    # We'll keep track of the current layer index (as a mutable list)
    current_index = [0]

    # A list to hold definitions for each layer.
    # Each entry is a dictionary with keys: name, polymer, thickness_value, thickness_unit, C0.
    layer_definitions = []
    def initialize_layer_definitions(n):
        nonlocal layer_definitions
        layer_definitions = []
        for i in range(n):
            layer_definitions.append({
                "name": f"P{i+1}",
                "polymer": default_polymer,
                "thickness_value": default_thickness_value,
                "thickness_unit": default_thickness_unit,
                "C0": default_c0
            })
    initialize_layer_definitions(num_layers_slider.value)

    # UI components for the current layer
    layer_name_input = widgets.Text(
        value=layer_definitions[0]["name"],
        description="Layer Name:",
        layout=widgets.Layout(width='50%')
    )
    polymer_dropdown = create_polymer_dropdown(default_value=default_polymer)
    thickness_value_input = widgets.FloatText(
         value=default_thickness_value,
         description='Thickness:',
         layout=widgets.Layout(width='30%')
    )
    thickness_unit_dropdown = widgets.Dropdown(
         options=["nm", "µm", "mm", "cm"],
         value=default_thickness_unit,
         description='Unit:',
         layout=widgets.Layout(width='20%')
    )
    c0_input = widgets.FloatText(
        value=default_c0,
        description="C0:",
        layout=widgets.Layout(width='30%')
    )

    # Navigation buttons: Previous and Next (using arrow icons)
    prev_button = widgets.Button(description="", icon="arrow-left")
    next_button = widgets.Button(description="", icon="arrow-right")

    # Button to instantiate all layers
    instantiate_all_button = widgets.Button(
         description="Instantiate All Layers",
         button_style="success",
         tooltip="Click to instantiate all layers"
    )

    # Output area for messages/results
    output = widgets.Output()

    # Function to load the current layer's definition into the UI fields
    def load_current_layer():
        idx = current_index[0]
        nav_label.value = f"Layer {idx+1} of {num_layers_slider.value}"
        data = layer_definitions[idx]
        layer_name_input.value = data["name"]
        polymer_dropdown.value = data["polymer"]
        thickness_value_input.value = data["thickness_value"]
        thickness_unit_dropdown.value = data["thickness_unit"]
        c0_input.value = data["C0"]

    # Function to save current UI field values into the current layer's definition
    def save_current_layer():
        idx = current_index[0]
        layer_definitions[idx]["name"] = layer_name_input.value.strip() or f"P{idx+1}"
        layer_definitions[idx]["polymer"] = polymer_dropdown.value
        layer_definitions[idx]["thickness_value"] = thickness_value_input.value
        layer_definitions[idx]["thickness_unit"] = thickness_unit_dropdown.value
        layer_definitions[idx]["C0"] = c0_input.value

    # Handlers for Previous and Next navigation buttons
    def on_prev(b):
        save_current_layer()
        if current_index[0] > 0:
            current_index[0] -= 1
            load_current_layer()

    def on_next(b):
        save_current_layer()
        if current_index[0] < num_layers_slider.value - 1:
            current_index[0] += 1
            load_current_layer()

    prev_button.on_click(on_prev)
    next_button.on_click(on_next)

    # Handler for number-of-layers slider change
    def on_num_layers_change(change):
        nonlocal layer_definitions
        if change['name'] == 'value':
            save_current_layer()
            new_n = change['new']
            current_n = len(layer_definitions)
            if new_n > current_n:
                # Append new layer definitions with defaults
                for i in range(current_n, new_n):
                    layer_definitions.append({
                        "name": f"P{i+1}",
                        "polymer": default_polymer,
                        "thickness_value": default_thickness_value,
                        "thickness_unit": default_thickness_unit,
                        "C0": default_c0
                    })
            elif new_n < current_n:
                layer_definitions = layer_definitions[:new_n]
                if current_index[0] >= new_n:
                    current_index[0] = new_n - 1
            load_current_layer()

    num_layers_slider.observe(on_num_layers_change, names='value')

    # Create a new text field for the assembly key (default "assembly")
    assembly_name_input = widgets.Text(
         value="multilayer1",
         description="Assembly Name:",
         layout=widgets.Layout(width='20%')
    )

    # Handler for the Instantiate All Layers button.
    def instantiate_all(b):
        save_current_layer()
        with output:
            output.clear_output()
            # Clear the global dictionary first.
            builtins.mylayers = {}
            # (We no longer clear a global mymaterial.)
            # Loop over all layer definitions and instantiate each layer.
            for i, layer_def in enumerate(layer_definitions):
                # Use the layer name (if empty, assign a default name)
                instance_name = layer_def["name"].strip() or f"P{i+1}"
                # Resolve the polymer class using the material surrogate
                PolymerClass = material(layer_def["polymer"])
                # Instantiate the layer; here D is fixed to 1e-12 as a demo.
                # Pass l and C0 as tuples; for C0, the unit is "a.u."
                layer_instance = PolymerClass(
                    l=(layer_def["thickness_value"], layer_def["thickness_unit"]),
                    D=1e-12,
                    C0=(layer_def["C0"], "a.u.")
                )
                builtins.mylayers[instance_name] = layer_instance
            # Instead of storing the sum in a single variable, store it in mymaterials using the assembly key.
            assembly_key = assembly_name_input.value.strip() or "assembly"
            # Here, layer.sum is assumed to be a class method of layer.
            builtins.mymaterials[assembly_key] = layer.sum(builtins.mylayers)
            print("Instantiated Layers:")
            for name, inst in builtins.mylayers.items():
                print(f"  {name}: {inst}")
            print(f"\nAssembly '{assembly_key}' has been created and stored in 'mymaterials'.")
    instantiate_all_button.on_click(instantiate_all)

    if _preheatedGUI_:
        instantiate_all(None) # we instantiate manually

    # Arrange the UI components into a layout.
    navigation_box = widgets.HBox([prev_button, nav_label, next_button])
    fields_box = widgets.VBox([
        layer_name_input,
        polymer_dropdown,
        widgets.HBox([thickness_value_input, thickness_unit_dropdown]),
        c0_input
    ])
    ui = widgets.VBox([
         num_layers_slider,
         navigation_box,
         fields_box,
         # Add the assembly name input above the instantiate button:
         assembly_name_input,
         instantiate_all_button,
         output
    ])

    # Load the first layer's data into the UI.
    load_current_layer()
    return ui



# %% layerLink
class layerLink:
    """
    A sparse representation of properties (`D`, `k`, `C0`) used in `layer` instances.

    This class allows storing and manipulating selected values of a property (`D`, `k`, or `C0`)
    while keeping a sparse structure. It enables seamless interaction with `layer` objects
    by overriding values dynamically and ensuring efficient memory usage.

    The primary use case is to fit and control property values externally while keeping
    the `layer` representation internally consistent.

    Attributes
    ----------
    property : str
        The name of the property linked (`"D"`, `"k"`, or `"C0"`).
    indices : np.ndarray
        A NumPy array storing the indices of explicitly defined values.
    values : np.ndarray
        A NumPy array storing the corresponding values at `indices`.
    length : int
        The total length of the sparse vector, ensuring coverage of all indices.
    replacement : str, optional
        Defines how missing values are handled:
        - `"repeat"`: Propagates the last known value beyond `length`.
        - `"periodic"`: Cycles through known values beyond `length`.
        - Default: No automatic replacement within `length`.

    Methods
    -------
    set(index, value)
        Sets values at specific indices. If `None` or `np.nan` is provided, the index is removed.
    get(index=None)
        Retrieves values at the given indices. Returns `NaN` for missing values.
    getandreplace(indices, altvalues)
        Similar to `get()`, but replaces `NaN` values with corresponding values from `altvalues`.
    getfull(altvalues)
        Returns the full vector using `getandreplace(None, altvalues)`.
    lengthextension()
        Ensures `length` covers all stored indices (`max(indices) + 1`).
    rename(new_property_name)
        Renames the `property` associated with this `layerLink`.
    nzcount()
        Returns the number of explicitly stored (nonzero) elements.
    __getitem__(index)
        Allows retrieval using `D_link[index]`, equivalent to `get(index)`.
    __setitem__(index, value)
        Allows assignment using `D_link[index] = value`, equivalent to `set(index, value)`.
    __add__(other)
        Concatenates two `layerLink` instances with the same property.
    __mul__(n)
        Repeats the `layerLink` instance `n` times, shifting indices accordingly.

    Examples
    --------
    Create a `layerLink` for `D` and manipulate its values:

    ```python
    D_link = layerLink("D")
    D_link.set([0, 2], [1e-14, 3e-14])
    print(D_link.get())  # Expected: array([1e-14, nan, 3e-14])

    D_link[1] = 2e-14
    print(D_link.get())  # Expected: array([1e-14, 2e-14, 3e-14])
    ```

    Concatenating two `layerLink` instances:

    ```python
    A = layerLink("D")
    A.set([0, 2], [1e-14, 3e-14])

    B = layerLink("D")
    B.set([1, 3], [2e-14, 4e-14])

    C = A + B  # Concatenates while shifting indices
    print(C.get())  # Expected: array([1e-14, 3e-14, nan, nan, 2e-14, 4e-14])
    ```

    Handling missing values with `getandreplace()`:

    ```python
    alt_values = np.array([5e-14, 6e-14, 7e-14, 8e-14])
    print(D_link.getandreplace([0, 1, 2, 3], alt_values))
    # Expected: array([1e-14, 2e-14, 3e-14, 8e-14])  # Fills NaNs from alt_values
    ```

    Ensuring correct behavior for `*`:

    ```python
    B = A * 3  # Repeats A three times
    print(B.indices)  # Expected: [0, 2, 4, 6, 8, 10]
    print(B.values)   # Expected: [1e-14, 3e-14, 1e-14, 3e-14, 1e-14, 3e-14]
    print(B.length)   # Expected: 3 * A.length
    ```


    Other Examples:
    ----------------

    ### **Creating a Link**
    D_link = layerLink("D", indices=[1, 3], values=[5e-14, 7e-14], length=4)
    print(D_link)  # <Link for D: 2 of 4 replacement values>

    ### **Retrieving Values**
    print(D_link.get())       # Full vector with None in unspecified indices
    print(D_link.get(1))      # Returns 5e-14
    print(D_link.get([0,2]))  # Returns [None, None]

    ### **Setting Values**
    D_link.set(2, 6e-14)
    print(D_link.get())  # Now index 2 is replaced

    ### **Resetting with a Prototype**
    prototype = [None, 5e-14, None, 7e-14, 8e-14]
    D_link.reset(prototype)
    print(D_link.get())  # Now follows the new structure

    ### **Getting and Setting Values with []**
    D_link = layerLink("D", indices=[1, 3, 5], values=[5e-14, 7e-14, 6e-14], length=10)
    print(D_link[3])      # ✅ Returns 7e-14
    print(D_link[:5])     # ✅ Returns first 5 elements (with NaNs where undefined)
    print(D_link[[1, 3]]) # ✅ Returns [5e-14, 7e-14]
    D_link[2] = 9e-14     # ✅ Sets D[2] to 9e-14
    D_link[0:4:2] = [1e-14, 2e-14]  # ✅ Sets D[0] = 1e-14, D[2] = 2e-14
    print(len(D_link))    # ✅ Returns 10 (full vector length)

    ###**Practical Syntaxes**
    D_link = layerLink("D")
    D_link[2] = 3e-14  # ✅ single value
    D_link[0] = 1e-14
    print(D_link.get())
    print(D_link[1])
    print(repr(D_link))
    D_link[:4] = 1e-16  # ✅ Fills indices 0,1,2,3 with 1e-16
    print(D_link.get())  # ✅ Outputs: [1e-16, 1e-16, 1e-16, 1e-16, nan, 1e-14]
    D_link[[1,2]] = None  # ✅ Fills indices 0,1,2,3 with 1e-16
    print(D_link.get())  # ✅ Outputs: [1e-16, 1e-16, 1e-16, 1e-16, nan, 1e-14]
    D_link[[0]] = 1e-10
    print(D_link.get())

    ###**How it works inside layer: a short simulation**
    # layerLink created by user
    duser = layerLink()
    duser.getfull([1e-15,2e-15,3e-15])
    duser[0] = 1e-10
    duser.getfull([1e-15,2e-15,3e-15])
    duser[1]=1e-9
    duser.getfull([1e-15,2e-15,3e-15])
    # layerLink used internally
    dalias=duser
    dalias[1]=2e-11
    duser.getfull([1e-15,2e-15,3e-15,4e-15])
    dalias[1]=2.1e-11
    duser.getfull([1e-15,2e-15,3e-15,4e-15])

    ###**Combining layerLinks instances**
    A = layerLink("D")
    A.set([0, 2], [1e-11, 3e-11])  # length=3
    B = layerLink("D")
    B.set([1, 3], [2e-14, 4e-12])  # length=4
    C = A + B
    print(C.indices)  # Expected: [0, 2, 4, 6]
    print(C.values)   # Expected: [1.e-11 3.e-11 2.e-14 4.e-12]
    print(C.length)   # Expected: 3 + 4 = 7


    TEST CASES:
    -----------

    print("🔹 Test 1: Initialize empty layerLink")
    D_link = layerLink("D")
    print(D_link.get())  # Expected: array([]) or array([nan, nan, nan]) if length is pre-set
    print(repr(D_link))  # Expected: No indices set

    print("\n🔹 Test 2: Assigning values at specific indices")
    D_link[2] = 3e-14
    D_link[0] = 1e-14
    print(D_link.get())  # Expected: array([1.e-14, nan, 3.e-14])
    print(D_link[1])     # Expected: nan

    print("\n🔹 Test 3: Assign multiple values at once")
    D_link[[1, 4]] = [2e-14, 5e-14]
    print(D_link.get())  # Expected: array([1.e-14, 2.e-14, 3.e-14, nan, 5.e-14])

    print("\n🔹 Test 4: Remove a single index")
    D_link[1] = None
    print(D_link.get())  # Expected: array([1.e-14, nan, 3.e-14, nan, 5.e-14])

    print("\n🔹 Test 5: Remove multiple indices at once")
    D_link[[0, 2]] = None
    print(D_link.get())  # Expected: array([nan, nan, nan, nan, 5.e-14])

    print("\n🔹 Test 6: Removing indices using a slice")
    D_link[3:5] = None
    print(D_link.get())  # Expected: array([nan, nan, nan, nan, nan])

    print("\n🔹 Test 7: Assign new values after removals")
    D_link[1] = 7e-14
    D_link[3] = 8e-14
    print(D_link.get())  # Expected: array([nan, 7.e-14, nan, 8.e-14, nan])

    print("\n🔹 Test 8: Check periodic replacement")
    D_link = layerLink("D", replacement="periodic")
    D_link[2] = 3e-14
    D_link[0] = 1e-14
    print(D_link[5])  # Expected: 1e-14 (since 5 mod 2 = 0)

    print("\n🔹 Test 9: Check repeat replacement")
    D_link = layerLink("D", replacement="repeat")
    D_link[2] = 3e-14
    D_link[0] = 1e-14
    print(D_link.get())  # Expected: array([1.e-14, nan, 3.e-14])
    print(D_link[3])     # Expected: 3e-14 (repeat last known value)

    print("\n🔹 Test 10: Resetting with a prototype")
    D_link.reset([None, 5e-14, None, 7e-14])
    print(D_link.get())  # Expected: array([nan, 5.e-14, nan, 7.e-14])

    print("\n🔹 Test 11: Edge case - Assigning nan explicitly")
    D_link[1] = np.nan
    print(D_link.get())  # Expected: array([nan, nan, nan, 7.e-14])

    print("\n🔹 Test 12: Assigning a range with a scalar value (broadcasting)")
    D_link[0:3] = 9e-14
    print(D_link.get())  # Expected: array([9.e-14, 9.e-14, 9.e-14, 7.e-14])

    print("\n🔹 Test 13: Assigning a slice with a list of values")
    D_link[1:4] = [6e-14, 5e-14, 4e-14]
    print(D_link.get())  # Expected: array([9.e-14, 6.e-14, 5.e-14, 4.e-14])

    print("\n🔹 Test 14: Length updates correctly after removals")
    D_link[[1, 2]] = None
    print(len(D_link))   # Expected: 4 (since max index is 3)

    print("\n🔹 Test 15: Setting index beyond length auto-extends")
    D_link[6] = 2e-14
    print(len(D_link))   # Expected: 7 (since max index is 6)
    print(D_link.get())  # Expected: array([9.e-14, nan, nan, 4.e-14, nan, nan, 2.e-14])

    """

    def __init__(self, property="D", indices=None, values=None, length=None,
                 replacement="repeat", dtype=np.float64, maxlength=None):
        """constructs a link"""
        self.property = property  # "D", "k", or "C0"
        self.replacement = replacement
        self.dtype = dtype
        self._maxlength = maxlength
        if isinstance(indices,(int,float)): indices = [indices]
        if isinstance(values,(int,float)): values = [values]

        if indices is None or values is None:
            self.indices = np.array([], dtype=int)
            self.values = np.array([], dtype=dtype)
        else:
            self.indices = np.array(indices, dtype=int)
            self.values = np.array(values, dtype=dtype)

        self.length = length if length is not None else (self.indices.max() + 1 if self.indices.size > 0 else 0)
        self._validate()

    def _validate(self):
        """Ensures consistency between indices and values."""
        if len(self.indices) != len(self.values):
            raise ValueError("indices and values must have the same length.")
        if self.indices.size > 0 and self.length < self.indices.max() + 1:
            raise ValueError("length must be at least max(indices) + 1.")

    def reset(self, prototypevalues):
        """
        Resets the link instance based on the prototype values.

        - Stores only non-None values.
        - Updates `indices`, `values`, and `length` accordingly.
        """
        self.indices = np.array([i for i, v in enumerate(prototypevalues) if v is not None], dtype=int)
        self.values = np.array([v for v in prototypevalues if v is not None], dtype=self.dtype)
        self.length = len(prototypevalues)  # Update the total length

    def get(self, index=None):
        """
        Retrieves values based on index or returns the full vector.

        Rules:
        - If `index=None`, returns the full vector with overridden values (no replacement applied).
        - If `index` is a scalar, returns the corresponding value, applying replacement rules if needed.
        - If `index` is an array, returns an array of the requested indices, applying replacement rules.

        Returns:
        - NumPy array with requested values.
        """
        if index is None:
            # Return the full vector WITHOUT applying any replacement
            full_vector = np.full(self.length, np.nan, dtype=self.dtype)
            full_vector[self.indices] = self.values  # Set known values
            return full_vector

        if np.isscalar(index):
            return self._get_single(index)

        # Ensure index is an array
        index = np.array(index, dtype=int)
        return np.array([self._get_single(i) for i in index], dtype=self.dtype)

    def _get_single(self, i):
        """Retrieves the value for a single index, applying rules if necessary."""
        if i in self.indices:
            return self.values[np.where(self.indices == i)[0][0]]

        if i >= self.length:  # Apply replacement *only* for indices beyond length
            if self.replacement == "periodic":
                return self.values[i % len(self.values)]
            elif self.replacement == "repeat":
                return self._get_single(self.length - 1)  # Repeat last known value

        return np.nan  # Default case for undefined in-bounds indices


    def set(self, index, value):
        """
        Sets values at specific indices.

        - If `index=None`, resets the link with `value`.
        - If `index` is a scalar, updates or inserts the value.
        - If `index` is an array, updates corresponding values.
        - If `value` is `None` or `np.nan`, removes the corresponding index.
        """
        if index is None:
            self.reset(value)
            return

        index = np.array(index, dtype=int)
        value = np.array(value, dtype=self.dtype)

        # check against _maxlength if defined
        if self._maxlength is not None:
            if np.any(index>=self._maxlength):
                raise IndexError(f"index cannot exceeds the number of layers-1 {self._maxlength-1}")

        # Handle scalars properly
        if np.isscalar(index):
            index = np.array([index])
            value = np.array([value])

        # Detect None or NaN values and remove those indices
        mask = np.isnan(value) if value.dtype.kind == 'f' else np.array([v is None for v in value])
        if np.any(mask):
            self._remove_indices(index[mask])  # Remove these indices
            index, value = index[~mask], value[~mask]  # Keep only valid values

        if index.size > 0:  # If there are remaining valid values, store them
            for i, v in zip(index, value):
                if i in self.indices:
                    self.values[np.where(self.indices == i)[0][0]] = v
                else:
                    self.indices = np.append(self.indices, i)
                    self.values = np.append(self.values, v)

        # Update length to ensure it remains valid
        if self.indices.size > 0:
            self.length = max(self.indices) + 1  # Adjust length based on max index
        else:
            self.length = 0  # Reset to 0 if empty

        self._validate()

    def _remove_indices(self, indices):
        """
        Removes indices from `self.indices` and `self.values` and updates length.
        """
        mask = np.isin(self.indices, indices, invert=True)
        self.indices = self.indices[mask]
        self.values = self.values[mask]

        # Update length after removal
        if self.indices.size > 0:
            self.length = max(self.indices) + 1  # Adjust length based on remaining max index
        else:
            self.length = 0  # Reset to 0 if no indices remain

    def reshape(self, new_length):
        """
        Reshapes the link instance to a new length.

        - If indices exceed new_length-1, they are removed with a warning.
        - If replacement operates beyond new_length-1, a warning is issued.
        """
        if new_length < self.length:
            invalid_indices = self.indices[self.indices >= new_length]
            if invalid_indices.size > 0:
                print(f"⚠️ Warning: Indices {invalid_indices.tolist()} are outside new length {new_length}. They will be removed.")
                mask = self.indices < new_length
                self.indices = self.indices[mask]
                self.values = self.values[mask]

        # Check if replacement would be applied beyond the new length
        if self.replacement == "repeat" and self.indices.size > 0 and self.length > new_length:
            print(f"⚠️ Warning: Repeat rule was defined for indices beyond {new_length-1}, but will not be used.")

        if self.replacement == "periodic" and self.indices.size > 0 and self.length > new_length:
            print(f"⚠️ Warning: Periodic rule was defined for indices beyond {new_length-1}, but will not be used.")

        self.length = new_length

    def __repr__(self):
        """Returns a detailed string representation."""
        txt = (f"Link(property='{self.property}', indices={self.indices.tolist()}, "
                f"values={self.values.tolist()}, length={self.length}, replacement='{self.replacement}')")
        print(txt)
        return(str(self))

    def __str__(self):
        """Returns a compact summary string."""
        return f"<{self.property}:{self.__class__.__name__}: {len(self.indices)}/{self.length}  values>"

    # Override `len()`
    def __len__(self):
        """Returns the length of the vector managed by the link object."""
        return self.length

    # Override `getitem` (support for indexing and slicing)
    def __getitem__(self, index):
        """
        Allows `D_link[index]` or `D_link[slice]` to retrieve values.

        - If `index` is an integer, returns a single value.
        - If `index` is a slice or list/array, returns a NumPy array of values.
        """
        if isinstance(index, slice):
            return self.get(np.arange(index.start or 0, index.stop or self.length, index.step or 1))
        return self.get(index)

    # Override `setitem` (support for indexing and slicing)
    def __setitem__(self, index, value):
        """
        Allows `D_link[index] = value` or `D_link[slice] = list/scalar`.

        - If `index` is an integer, updates or inserts a single value.
        - If `index` is a slice or list/array, updates multiple values.
        - If `value` is `None` or `np.nan`, removes the corresponding index.
        """
        if isinstance(index, slice):
            indices = np.arange(index.start or 0, index.stop or self.length, index.step or 1)

        elif isinstance(index, (list, np.ndarray)):  # Handle non-contiguous indices
            indices = np.array(index, dtype=int)

        elif np.isscalar(index):  # Single index assignment
            indices = np.array([index], dtype=int)

        else:
            raise TypeError(f"Unsupported index type: {type(index)}")

        if value is None or (isinstance(value, float) and np.isnan(value)):  # Remove these indices
            self._remove_indices(indices)
        else:
            values = np.full_like(indices, value, dtype=self.dtype) if np.isscalar(value) else np.array(value, dtype=self.dtype)
            if len(indices) != len(values):
                raise ValueError(f"Cannot assign {len(values)} values to {len(indices)} indices.")
            self.set(indices, values)

    def getandreplace(self, indices=None, altvalues=None):
        """
        Retrieves values for the given indices, replacing NaN values with corresponding values from altvalues.

        - If `indices` is None or empty, it defaults to `[0, 1, ..., self.length - 1]`
        - altvalues should be a NumPy array with the same dtype as self.values.
        - altvalues **can be longer than** self.length, but **cannot be shorter than the highest requested index**.
        - If an index is undefined (`NaN` in get()), it is replaced with altvalues[index].

        Parameters:
        ----------
        indices : list or np.ndarray (default: None)
            The indices to retrieve values for. If None, defaults to full range `[0, ..., self.length - 1]`.
        altvalues : list or np.ndarray
            Alternative values to use where `get()` returns `NaN`.

        Returns:
        -------
        np.ndarray
            A NumPy array of values, with NaNs replaced by altvalues.
        """
        if indices is None or len(indices) == 0:
            indices = np.arange(self.length)  # Default to full range

        indices = np.array(indices, dtype=int)
        altvalues = np.array(altvalues, dtype=self.dtype)

        max_requested_index = indices.max() if indices.size > 0 else 0
        if max_requested_index >= altvalues.shape[0]:  # Ensure altvalues covers all requested indices
            raise ValueError(
                f"altvalues is too short! It has length {altvalues.shape[0]}, but requested index {max_requested_index}."
            )
        # Get original values
        original_values = self.get(indices)
        # Replace NaN values with corresponding values from altvalues
        mask_nan = np.isnan(original_values)
        original_values[mask_nan] = altvalues[indices[mask_nan]]
        return original_values


    def getfull(self, altvalues):
        """
        Retrieves the full vector using `getandreplace(None, altvalues)`.

        - If `length == 0`, returns `altvalues` as a NumPy array of the correct dtype.
        - Extends `self.length` to match `altvalues` if it's shorter.
        - Supports multidimensional `altvalues` by flattening it.

        Parameters:
        ----------
        altvalues : list or np.ndarray
            Alternative values to use where `get()` returns `NaN`.

        Returns:
        -------
        np.ndarray
            Full vector with NaNs replaced by altvalues.
        """
        # Convert altvalues to a NumPy array and flatten if needed
        altvalues = np.array(altvalues, dtype=self.dtype).flatten()

        # If self has no length, return altvalues directly
        if self.length == 0:
            return altvalues

        # Extend self.length to match altvalues if needed
        if self.length < altvalues.shape[0]:
            self.length = altvalues.shape[0]

        return self.getandreplace(None, altvalues)

    @property
    def nzlength(self):
        """
        Returns the number of stored nonzero elements (i.e., indices with values).
        """
        return len(self.indices)

    def lengthextension(self):
        """
        Ensures that the length of the layerLink instance is at least `max(indices) + 1`.

        - If there are no indices, the length remains unchanged.
        - If `length` is already sufficient, nothing happens.
        - Otherwise, it extends `length` to `max(indices) + 1`.
        """
        if self.indices.size > 0:  # Only extend if there are indices
            self.length = max(self.length, max(self.indices) + 1)

    def rename(self, new_property_name):
        """
        Renames the property associated with this link.

        Parameters:
        ----------
        new_property_name : str
            The new property name.

        Raises:
        -------
        TypeError:
            If `new_property_name` is not a string.
        """
        if not isinstance(new_property_name, str):
            raise TypeError(f"Property name must be a string, got {type(new_property_name).__name__}.")
        self.property = new_property_name


    def __add__(self, other):
        """
        Concatenates two layerLink instances.

        - Only allowed if both instances have the same property.
        - Calls `lengthextension()` on both instances before summing lengths.
        - Shifts `other`'s indices by `self.length` to maintain sparsity.
        - Concatenates values and indices.

        Returns:
        -------
        layerLink
            A new concatenated layerLink instance.
        """
        if not isinstance(other, layerLink):
            raise TypeError(f"Cannot concatenate {type(self).__name__} with {type(other).__name__}")

        if self.property != other.property:
            raise ValueError(f"Cannot concatenate: properties do not match ('{self.property}' vs. '{other.property}')")

        # Ensure lengths are properly extended before computing new length
        self.lengthextension()
        other.lengthextension()

        # Create a new instance for the result
        result = layerLink(self.property)

        # Copy self's values
        result.indices = np.array(self.indices, dtype=int)
        result.values = np.array(self.values, dtype=self.dtype)

        # Adjust other’s indices and add them
        shifted_other_indices = np.array(other.indices) + self.length
        result.indices = np.concatenate([result.indices, shifted_other_indices])
        result.values = np.concatenate([result.values, np.array(other.values, dtype=self.dtype)])

        # ✅ Correct length calculation: Sum of the two lengths (assuming lengths are extended)
        result.length = self.length + other.length

        return result


    def __mul__(self, n):
        """
        Repeats the layerLink instance `n` times.

        - Uses `+` to concatenate multiple copies with shifted indices.
        - Each repetition gets indices shifted by `self.length * i`.

        Returns:
        -------
        layerLink
            A new layerLink instance with repeated data.
        """
        if not isinstance(n, int) or n <= 0:
            raise ValueError("Multiplication factor must be a positive integer")

        result = layerLink(self.property)
        for i in range(n):
            shifted_instance = layerLink(self.property)
            shifted_instance.indices = np.array(self.indices) + i * self.length
            shifted_instance.values = np.array(self.values, dtype=self.dtype)
            shifted_instance.length = self.length
            result += shifted_instance  # Use `+` to merge each repetition

        return result

# %% Core class: layer
# default values (usable in layer object methods)
# these default values can be moved in a configuration file


# Main class definition
# =======================
class layer:
    """
    ------------------------------------------------------------------------------
    **Core Functionality**
    ------------------------------------------------------------------------------
    This class models layers in food packaging, handling mass transfer, partitioning,
    and meshing for finite-volume simulations using a modified Patankar method.
    Layers can be assembled into multilayers via the `+` operator and support
    dynamic property linkage using `layerLink`.

    ------------------------------------------------------------------------------
    **Key Properties**
    ------------------------------------------------------------------------------
    - `l`: Thickness of the layer (m)
    - `D`: Diffusion coefficient (m²/s)
    - `k`: Partition coefficient (dimensionless)
    - `C0`: Initial concentration (arbitrary units)
    - `rho`: Density (kg/m³)
    - `T`: Contact temperature (°C)
    - `substance`: Migrant/substance modeled for diffusion
    - `medium`: The food medium in contact with the layer
    - `Dmodel`, `kmodel`: Callable models for diffusion and partitioning

    ------------------------------------------------------------------------------
    **Methods**
    ------------------------------------------------------------------------------
    - `__add__(self, other)`: Combines two layers into a multilayer structure.
    - `__mul__(self, n)`: Duplicates a layer `n` times to create a multilayer.
    - `__getitem__(self, i)`: Retrieves a sublayer from a multilayer.
    - `__setitem__(self, i, other)`: Replaces sublayers in a multilayer structure.
    - `mesh(self)`: Generates a numerical mesh for finite-volume simulations.
    - `struct(self)`: Returns a dictionary representation of the layer properties.
    - `resolvename(param_value, param_key, **unresolved)`: Resolves synonyms for parameter names.
    - `help(cls)`: Displays a dynamically formatted summary of input parameters.

    ------------------------------------------------------------------------------
    **Integration with SFPPy Modules**
    ------------------------------------------------------------------------------
    - Works with `migration.py` for mass transfer simulations.
    - Interfaces with `food.py` to define food-contact conditions.
    - Uses `property.py` for predicting diffusion (`D`) and partitioning (`k`).
    - Connects with `geometry.py` for 3D packaging simulations.

    ------------------------------------------------------------------------------
    **Usage Example**
    ------------------------------------------------------------------------------
    ```python
    from patankar.layer import LDPE, PP, layerLink

    # Define a polymer layer with default properties
    A = LDPE(l=50e-6, D=1e-14)

    # Create a multilayer structure
    B = PP(l=200e-6, D=1e-15)
    multilayer = A + B

    # Assign dynamic property linkage
    k_link = layerLink("k", indices=[1], values=[10])  # Assign partition coefficient to the second layer
    multilayer.klink = k_link

    # Simulate migration
    from patankar.migration import senspatankar
    from patankar.food import ethanol
    medium = ethanol()
    solution = senspatankar(multilayer, medium)
    solution.plotCF()
    ```

    ------------------------------------------------------------------------------
    **Notes**
    ------------------------------------------------------------------------------
    - This class supports dynamic property inheritance, meaning `D` and `k` can be computed
      based on the substance defined in `substance` and `medium`.
    - The `layerLink` mechanism allows parameter adjustments without modifying the core object.
    - The modified finite-volume meshing ensures **accurate steady-state and transient** behavior.

    """

    # -----------------------------------------------------------------------------
    # Class attributes that can be overidden in instances.
    # Their default values are set in classes and overriden with similar
    # instance properties with @property.setter.
    # These values cannot be set during construction, but only after instantiation.
    # -----------------------------------------------------------------------------
    # These properties are essential for model predictions, they cannot be customized
    # beyond the rules accepted by the model predictors (they are not metadata)
    _physicalstate = "solid"        # solid (default), liquid, gas, porous
    _chemicalclass = "polymer"      # polymer (default), other
    _chemicalsubstance = None       # None (default), monomer for polymers
    _polarityindex = 0.0            # polarity index (roughly: 0=hexane, 10=water)

    # Low-level prediction properties (these properties are common with patankar.food)
    _lowLevelPredictionPropertyList = ["physicalstate","chemicalclass",
                                       "chemicalsubstance","polarityindex","ispolymer","issolid"]

    # --------------------------------------------------------------------
    # PRIVATE PROPERTIES (cannot be changed by the user)
    # __ read only attributes
    #  _ private attributes (not public)
    # --------------------------------------------------------------------
    __description = "LAYER object"                # description
    __version = 1.0                               # version
    __contact = "olivier.vitrac@agroparistech.fr" # contact person
    _printformat = "%0.4g"   # format to display D, k, l values


    # Synonyms dictionary: Maps alternative names to the actual parameter
    # these synonyms can be used during construction
    _synonyms = {
        "substance": {"migrant", "compound", "chemical","molecule","solute"},
        "medium": {"food","simulant","fluid","liquid","contactmedium"},
        "C0": {"CP0", "Cp0"},
        "l": {"lp", "lP"},
        "D": {"Dp", "DP"},
        "k": {"kp", "kP"},
        "T": {"temp","Temp","temperature","Temperature",
              "contacttemperature","ContactTemperature","contactTemperature"}
    }
    # Default values for parameters (note that Td cannot be changed by the end-user)
    _defaults = {
        "l": 5e-5,   # Thickness (m)
        "D": 1e-14,  # Diffusion coefficient (m^2/s)
        "k": 1.0,      # Henri-like coefficient (dimensionless)
        "C0": 1000,  # Initial concentration (arbitrary units)
        "rho": 1000, # Default density (kg/m³)
        "T": 40.0,     # Default temperature (°C)
        "Td": 25.0,    # Reference temperature for densities (°C)
        # Units (do not change)
        "lunit": "m",
        "Dunit": "m**2/s",
        "kunit": "a.u.",  # NoUnits
        "Cunit": "a.u.",  # NoUnits
        "rhounit": "kg/m**3",
        "Tunit": "degC",  # Temperatures are indicated in °C instead of K (to reduce end-user mistakes)
        # Layer properties
        "layername": "my layer",
        "layertype": "unknown type",
        "layermaterial": "unknown material",
        "layercode": "N/A",
        # Mesh parameters
        "nmeshmin": 20,
        "nmesh": 600,
        # Substance
        "substance": None,
        "simulant": None,
        # Other parameters
        "verbose": None,
        "verbosity": 2
    }

    # List units
    _parametersWithUnits = {
        "l": "m",
        "D": "m**2/s",
        "k": "a.u.",
        "C": "a.u.",
        "rhp": "kg/m**3",
        "T": "degC",
        }

    # Brief descriptions for each parameter
    _descriptionInputs = {
        "l": "Thickness of the layer (m)",
        "D": "Diffusion coefficient (m²/s)",
        "k": "Henri-like coefficient (dimensionless)",
        "C0": "Initial concentration (arbitrary units)",
        "rho": "Density of the material (kg/m³)",
        "T": "Layer temperature (°C)",
        "Td": "Reference temperature for densities (°C)",
        "lunit": "Unit of thickness (default: m)",
        "Dunit": "Unit of diffusion coefficient (default: m²/s)",
        "kunit": "Unit of Henri-like coefficient (default: a.u.)",
        "Cunit": "Unit of initial concentration (default: a.u.)",
        "rhounit": "Unit of density (default: kg/m³)",
        "Tunit": "Unit of temperature (default: degC)",
        "layername": "Name of the layer",
        "layertype": "Type of layer (e.g., polymer, ink, air)",
        "layermaterial": "Material composition of the layer",
        "layercode": "Identification code for the layer",
        "nmeshmin": "Minimum number of FV mesh elements for the layer",
        "nmesh": "Number of FV mesh elements for numerical computation",
        "verbose": "Verbose mode (None or boolean)",
        "verbosity": "Level of verbosity for debug messages (integer)"
    }

    # --------------------------------------------------------------------
    # CONSTRUCTOR OF INSTANCE PROPERTIES
    # None = missing numeric value (managed by default)
    # --------------------------------------------------------------------
    def __init__(self,
                 l=None, D=None, k=None, C0=None, rho=None, T=None,
                 lunit=None, Dunit=None, kunit=None, Cunit=None, rhounit=None, Tunit=None,
                 layername=None,layertype=None,layermaterial=None,layercode=None,
                 substance = None, medium = None,
                 # Dmodel = None, kmodel = None, they are defined via migrant (future overrides)
                 nmesh=None, nmeshmin=None, # simulation parametes
                 # link properties (for fitting and linking properties across simulations)
                 Dlink=None, klink=None, C0link=None, Tlink=None, llink=None,
                 verbose=None, verbosity=2,**unresolved):
        """

        Parameters
        ----------

        layername : TYPE, optional, string
                    DESCRIPTION. Layer Name. The default is "my layer".
        layertype : TYPE, optional, string
                    DESCRIPTION. Layer Type. The default is "unknown type".
        layermaterial : TYPE, optional, string
                        DESCRIPTION. Material identification . The default is "unknown material".
        PHYSICAL QUANTITIES
        l : TYPE, optional, scalar or tupple (value,"unit")
            DESCRIPTION. Thickness. The default is 50e-6 (m).
        D : TYPE, optional, scalar or tupple (value,"unit")
            DESCRIPTION. Diffusivity. The default is 1e-14 (m^2/s).
        k : TYPE, optional, scalar or tupple (value,"unit")
            DESCRIPTION. Henry-like coefficient. The default is 1 (a.u.).
        C0 : TYPE, optional, scalar or tupple (value,"unit")
            DESCRIPTION. Initial concentration. The default is 1000 (a.u.).
        PHYSICAL UNITS
        lunit : TYPE, optional, string
                DESCRIPTION. Length units. The default unit is "m.
        Dunit : TYPE, optional, string
                DESCRIPTION. Diffusivity units. The default unit is 1e-14 "m^2/s".
        kunit : TYPE, optional, string
                DESCRIPTION. Henry-like coefficient. The default unit is "a.u.".
        Cunit : TYPE, optional, string
                DESCRIPTION. Initial concentration. The default unit is "a.u.".
        Returns
        -------
        a monolayer object which can be assembled into a multilayer structure

        """
        # resolve alternative names used by end-users
        substance = layer.resolvename(substance,"substance",**unresolved)
        medium = layer.resolvename(medium, "medium", **unresolved)
        C0 = layer.resolvename(C0,"C0",**unresolved)
        l = layer.resolvename(l,"l",**unresolved)
        D = layer.resolvename(D,"D",**unresolved)
        k = layer.resolvename(k,"k",**unresolved)
        T = layer.resolvename(T,"T",**unresolved)

        # Assign defaults only if values are not provided
        l = l if l is not None else layer._defaults["l"]
        D = D if D is not None else layer._defaults["D"]
        k = k if k is not None else layer._defaults["k"]
        C0 = C0 if C0 is not None else layer._defaults["C0"]
        rho = rho if rho is not None else layer._defaults["rho"]
        T = T if T is not None else layer._defaults["T"]
        lunit = lunit if lunit is not None else layer._defaults["lunit"]
        Dunit = Dunit if Dunit is not None else layer._defaults["Dunit"]
        kunit = kunit if kunit is not None else layer._defaults["kunit"]
        Cunit = Cunit if Cunit is not None else layer._defaults["Cunit"]
        rhounit = rhounit if rhounit is not None else layer._defaults["rhounit"]
        Tunit = Tunit if Tunit is not None else layer._defaults["Tunit"]
        nmesh = nmesh if nmesh is not None else layer._defaults["nmesh"]
        nmeshmin = nmeshmin if nmeshmin is not None else layer._defaults["nmeshmin"]
        verbose = verbose if verbose is not None else layer._defaults["verbose"]
        verbosity = verbosity if verbosity is not None else layer._defaults["verbosity"]

        # add user overrides
        nmesh = useroverride("nmesh",nmesh)
        nmeshmin = useroverride("nmeshin",nmeshmin)

        # Assign layer id properties
        layername = layername if layername is not None else layer._defaults["layername"]
        layertype = layertype if layertype is not None else layer._defaults["layertype"]
        layermaterial = layermaterial if layermaterial is not None else layer._defaults["layermaterial"]
        layercode = layercode if layercode is not None else layer._defaults["layercode"]

        # validate all physical paramaters with their units
        l,lunit = check_units(l,lunit,layer._defaults["lunit"])
        D,Dunit = check_units(D,Dunit,layer._defaults["Dunit"])
        k,kunit = check_units(k,kunit,layer._defaults["kunit"])
        C0,Cunit = check_units(C0,Cunit,layer._defaults["Cunit"])
        rho,rhounit = check_units(rho,rhounit,layer._defaults["rhounit"])
        T,Tunit = check_units(T,Tunit,layer._defaults["Tunit"])

        # set attributes: id and physical properties
        self._name = [layername]
        self._type = [layertype]
        self._material = [layermaterial]
        self._code = [layercode]
        self._nlayer = 1
        self._l = l[:1]
        self._D = D[:1]
        self._k = k[:1]
        self._C0 = C0[:1]
        self._rho = rho[:1]
        self._T = T
        self._lunit = lunit
        self._Dunit = Dunit
        self._kunit = kunit
        self._Cunit = Cunit
        self._rhounit = rhounit
        self._Tunit = Tunit
        self._nmesh = nmesh
        self._nmeshmin = nmeshmin

        # intialize links for X = D,k,C0,T,l (see documentation of layerLink)
        # A link enables the values of X to be defined and controlled outside the instance
        self._Dlink  = self._initialize_link(Dlink, "D")
        self._klink  = self._initialize_link(klink, "k")
        self._C0link = self._initialize_link(C0link, "C0")
        self._Tlink  = self._initialize_link(Tlink, "T")
        self._llink  = self._initialize_link(llink, "l")

        # set substance, medium and related D and k models
        if isinstance(substance,str):
            substance = migrant(substance)
        if substance is not None and not isinstance(substance,migrant):
            raise ValueError(f"subtance must be None a or a migrant not a {type(substance).__name__}")
        self._substance = substance
        if medium is not None:
            from patankar.food import foodlayer # local import only if needed
            if not isinstance(medium,foodlayer):
                raise ValueError(f"medium must be None or a foodlayer not a {type(medium).__name__}")
        self._medium = medium
        self._Dmodel = "default"  # do not use directly self._compute_Dmodel (force refresh)
        self._kmodel = "default"  # do not use directly self._compute_kmodel (force refresh)

        # set history for all layers merged with +
        self._layerclass_history = []
        self._ispolymer_history = []
        self._chemicalsubstance_history = []
        self._Tg_history = []
        self._porosity_history = []
        self._crystallinity_history = []

        # set verbosity attributes
        self.verbosity = 0 if verbosity is None else verbosity
        self.verbose = verbosity>0 if verbose is None else verbose

        # we initialize the acknowlegment process for future property propagation
        self._hasbeeninherited = {}


    # --------------------------------------------------------------------
    # Helper method: initializes and validates layerLink attributes
    # (Dlink, klink, C0link, Tlink, llink)
    # --------------------------------------------------------------------
    def _initialize_link(self, link, expected_property):
        """
        Initializes and validates a layerLink attribute.

        Parameters:
        ----------
        link : layerLink or None
            The `layerLink` instance to be assigned.
        expected_property : str
            The expected property name (e.g., "D", "k", "C0", "T").

        Returns:
        -------
        layerLink or None
            The validated `layerLink` instance or None.

        Raises:
        -------
        TypeError:
            If `link` is not a `layerLink` or `None`.
        ValueError:
            If `link.property` does not match `expected_property`.
        """
        if link is None:
            return None
        if isinstance(link, layerLink):
            if link.property == expected_property:
                return link
            raise ValueError(f'{expected_property}link.property should be "{expected_property}" not "{link.property}"')
        raise TypeError(f"{expected_property}link must be a layerLink not a {type(link).__name__}")


    # --------------------------------------------------------------------
    # Class method returning help() for the end user
    # --------------------------------------------------------------------
    @classmethod
    def help(cls):
        """
        Prints a dynamically formatted summary of all input parameters,
        adjusting column widths based on content and wrapping long descriptions.
        """

        # Column Headers
        headers = ["Parameter", "Default Value", "Has Synonyms?", "Description"]
        col_widths = [len(h) for h in headers]  # Start with header widths

        # Collect Data Rows
        rows = []
        for param, default in cls._defaults.items():
            has_synonyms = "✅ Yes" if param in cls._synonyms else "❌ No"
            description = cls._descriptionInputs.get(param, "No description available")

            # Update column widths dynamically
            col_widths[0] = max(col_widths[0], len(param))
            col_widths[1] = max(col_widths[1], len(str(default)))
            col_widths[2] = max(col_widths[2], len(has_synonyms))
            col_widths[3] = max(col_widths[3], len(description))

            rows.append([param, str(default), has_synonyms, description])

        # Function to wrap text for a given column width
        def wrap_text(text, width):
            return textwrap.fill(text, width)

        # Print Table with Adjusted Column Widths
        separator = "+-" + "-+-".join("-" * w for w in col_widths) + "-+"
        print("\n### **Accepted Parameters and Defaults**\n")
        print(separator)
        print("| " + " | ".join(h.ljust(col_widths[i]) for i, h in enumerate(headers)) + " |")
        print(separator)
        for row in rows:
            # Wrap text in the description column
            row[3] = wrap_text(row[3], col_widths[3])

            # Print row
            print("| " + " | ".join(row[i].ljust(col_widths[i]) for i in range(3)) + " | " + row[3])
        print(separator)

        # Synonyms Table
        print("\n### **Parameter Synonyms**\n")
        syn_headers = ["Parameter", "Synonyms"]
        syn_col_widths = [
            max(len("Parameter"), max(len(k) for k in cls._synonyms.keys())),  # Ensure it fits "Parameter"
            max(len("Synonyms"), max(len(", ".join(v)) for v in cls._synonyms.values()))  # Ensure it fits "Synonyms"
        ]
        syn_separator = "+-" + "-+-".join("-" * w for w in syn_col_widths) + "-+"
        print(syn_separator)
        print("| " + " | ".join(h.ljust(syn_col_widths[i]) for i, h in enumerate(syn_headers)) + " |")
        print(syn_separator)
        for param, synonyms in cls._synonyms.items():
            print(f"| {param.ljust(syn_col_widths[0])} | {', '.join(synonyms).ljust(syn_col_widths[1])} |")
        print(syn_separator)


    # --------------------------------------------------------------------
    # Class method to handle ambiguous definitions from end-user
    # --------------------------------------------------------------------
    @classmethod
    def resolvename(cls, param_value, param_key, **unresolved):
        """
        Resolves the correct parameter value using known synonyms.

        - If param_value is already set (not None), return it.
        - If a synonym exists in **unresolved, assign its value.
        - If multiple synonyms of the same parameter appear in **unresolved, raise an error.
        - Otherwise, return None.

        Parameters:
        - `param_name` (any): The original value (if provided).
        - `param_key` (str): The legitimate parameter name we are resolving.
        - `unresolved` (dict): The dictionary of unrecognized keyword arguments.

        Returns:
        - The resolved value or None if not found.
        """
        if param_value is not None:
            return param_value  # The parameter is explicitly defined, do not override
        if not unresolved:      # shortcut
            return None
        resolved_value = None
        found_keys = []
        # Check if param_key itself is present in unresolved
        if param_key in unresolved:
            found_keys.append(param_key)
            resolved_value = unresolved[param_key]
        # Check if any of its synonyms are in unresolved
        if param_key in cls._synonyms:
            for synonym in cls._synonyms[param_key]:
                if synonym in unresolved:
                    found_keys.append(synonym)
                    resolved_value = unresolved[synonym]
        # Raise error if multiple synonyms were found
        if len(found_keys) > 1:
            raise ValueError(
                f"Conflicting definitions: Multiple synonyms {found_keys} were provided for '{param_key}'."
            )
        return resolved_value


    # --------------------------------------------------------------------
    # overloading binary addition (note that the output is of type layer)
    # --------------------------------------------------------------------
    def __add__(self, other):
        """ C = A + B | overload + operator """
        if isinstance(other, layer):
            res = duplicate(self)
            res._nmeshmin = min(self._nmeshmin, other._nmeshmin)
            # Propagate substance
            if self._substance is None:
                res._substance = other._substance
            else:
                if isinstance(self._substance, migrant) and isinstance(other._substance, migrant):
                    if self._substance.M != other._substance.M:
                        print("Warning: the smallest substance is propagated everywhere")
                    res._substance = self._substance if self._substance.M <= other._substance.M else other._substance
                else:
                    res._substance = None
            # Concatenate general attributes
            for p in ["_name", "_type", "_material", "_code", "_nlayer"]:
                setattr(res, p, getattr(self, p) + getattr(other, p))
            # Concatenate numeric arrays
            for p in ["_l", "_D", "_k", "_C0", "_rho", "_T"]:
                setattr(res, p, np.concatenate((getattr(self, p), getattr(other, p))))
            # Handle history tracking
            res._layerclass_history = self.layerclass_history + other.layerclass_history
            res._ispolymer_history = self.ispolymer_history + other.ispolymer_history
            res._chemicalsubstance_history = self.chemicalsubstance_history + other.chemicalsubstance_history
            res._Tg_history = self.Tg_history + other.Tg_history
            res._porosity_history = self.porosity_history + other.porosity_history
            res._crystallinity_history = self.crystallinity_history + other.crystallinity_history
            # Manage layerLink attributes (Dlink, klink, C0link, Tlink, llink)
            property_map = {
                "Dlink": ("D", self.Dlink, other.Dlink),
                "klink": ("k", self.klink, other.klink),
                "C0link": ("C0", self.C0link, other.C0link),
                "Tlink": ("T", self.Tlink, other.Tlink),
                "llink": ("l", self.llink, other.llink),
            }
            for attr, (prop, self_link, other_link) in property_map.items():
                if (self_link is not None) and (other_link is not None):
                    # Case 1: Both have a link → Apply `+`
                    setattr(res, '_'+attr, self_link + other_link)
                elif self_link is not None:
                    # Case 2: Only self has a link → Use as-is
                    setattr(res, '_'+attr, self_link)
                elif other_link is not None:
                    # Case 3: Only other has a link → Shift indices and use
                    shifted_link = duplicate(other_link)
                    shifted_link.indices += len(getattr(self, prop))
                    setattr(res, '_'+attr, shifted_link)
                else:
                    # Case 4: Neither has a link → Result is None
                    setattr(res, '_'+attr, None)
            return res
        else:
            raise ValueError("Invalid layer object")

    # sum: + over a list/dict
    @classmethod
    def sum(cls, layers):
        """
        Sums all layer instances from a dictionary or list and returns a multilayer instance.

        Parameters:
          layers (dict or list): A dictionary (whose values are layer instances) or a list of layer instances.

        Returns:
          A single layer instance representing the stack (sum) of all input layers.

        Example:
          mymultilayer = layer.sum(mylayers)
          # where 'mylayers' could be a dict of layers or a list of layers.
        """
        # If layers is a dictionary, iterate over its values.
        if isinstance(layers, dict):
            layer_iter = layers.values()
        else:
            layer_iter = layers
        total = None
        for L in layer_iter:
            if not isinstance(L,layer):
                raise TypeError(f"all layers must of class layer not a {type(L).__name__}")
            if total is None:
                total = L
            else:
                total = total + L  # This uses the overloaded __add__ operator of the layer class.
        return total


    # --------------------------------------------------------------------
    # overloading binary multiplication (note that the output is of type layer)
    # --------------------------------------------------------------------
    def __mul__(self,ntimes):
        """ nA = A*n | overload * operator """
        if isinstance(ntimes, int) and ntimes>0:
            res = duplicate(self)
            if ntimes>1:
                for n in range(1,ntimes): res += self
            return res
        else: raise ValueError("multiplicator should be a strictly positive integer")


    # --------------------------------------------------------------------
    # len method
    # --------------------------------------------------------------------
    def __len__(self):
        """ length method """
        return self._nlayer

    # --------------------------------------------------------------------
    # object indexing (get,set) method
    # --------------------------------------------------------------------
    def __getitem__(self,i):
        """ get indexing method """
        res = duplicate(self)
        # check indices
        isscalar = isinstance(i,int)
        if isinstance(i,slice):
            if i.step==None: j = list(range(i.start,i.stop))
            else: j = list(range(i.start,i.stop,i.step))
            res._nlayer = len(j)
        if isinstance(i,int): res._nlayer = 1
        # pick indices for each property
        for p in ["_name","_type","_material","_l","_D","_k","_C0",
                  # Handle history tracking
                  "_layerclass_history","_ispolymer_history","_chemicalsubstance_history",
                  "_Tg_history","_porosity_history","_crystallinity_history"]:
            content = getattr(self,p)
            try:
                if isscalar: setattr(res,p,content[i:i+1])
                else: setattr(res,p,content[i])
            except IndexError as err:
                if self.verbosity>0 and self.verbose:
                    print("bad layer object indexing: ",err)
        return res

    def __setitem__(self,i,other):
        """ set indexing method """
        # check indices
        if isinstance(i,slice):
            if i.step==None: j = list(range(i.start,i.stop))
            else: j = list(range(i.start,i.stop,i.step))
        elif isinstance(i,int): j = [i]
        else:raise IndexError("invalid index")
        islayer = isinstance(other,layer)
        isempty = not islayer and isinstance(other,list) and len(other)<1
        if isempty:         # empty right hand side
            for p in ["_name","_type","_material","_l","_D","_k","_C0",
                      # Handle history tracking
                      "_layerclass_history","_ispolymer_history","_chemicalsubstance_history",
                      "_Tg_history","_porosity_history","_crystallinity_history"]:
                content = getattr(self,p)
                try:
                    newcontent = [content[k] for k in range(self._nlayer) if k not in j]
                except IndexError as err:
                    if self.verbosity>0 and self.verbose:
                        print("bad layer object indexing: ",err)
                if isinstance(content,np.ndarray) and not isinstance(newcontent,np.ndarray):
                    newcontent = np.array(newcontent)
                setattr(self,p,newcontent)
            self._nlayer = len(newcontent)
        elif islayer:        # islayer right hand side
            nk1 = len(j)
            nk2 = other._nlayer
            if nk1 != nk2:
                raise IndexError("the number of elements does not match the number of indices")
            for p in ["_name","_type","_material","_l","_D","_k","_C0"
                      # Handle history tracking
                      "_layerclass_history","_ispolymer_history","_chemicalsubstance_history",
                      "_Tg_history","_porosity_history","_crystallinity_history"]:
                content1 = getattr(self,p)
                content2 = getattr(other,p)
                for k in range(nk1):
                    try:
                        content1[j[k]] = content2[k]
                    except IndexError as err:
                        if self.verbosity>0 and self.verbose:
                            print("bad layer object indexing: ",err)
                setattr(self,p,content1)
        else:
            raise ValueError("only [] or layer object are accepted")


    # --------------------------------------------------------------------
    # Getter methods (show private/hidden properties and meta-properties)
    # --------------------------------------------------------------------
    # Return class or instance attributes
    @property
    def physicalstate(self): return self._physicalstate
    @property
    def chemicalclass(self): return self._chemicalclass
    @property
    def chemicalsubstance(self): return self._chemicalsubstance
    @property
    def polarityindex(self):
        # rescaled to match predictions - standard scale [0,10.2] - predicted scale [0,7.12]
        return self._polarityindex * migrant("water").polarityindex/10.2
    @property
    def ispolymer(self): return self.chemicalclass == "polymer"
    @property
    def issolid(self): return self.physicalstate == "solid"
    @property
    def layerclass_history(self):
        return self._layerclass_history if self._layerclass_history != [] else [self.layerclass]
    @property
    def ispolymer_history(self):
        return self._ispolymer_history if self._ispolymer_history != [] else [self.ispolymer]
    @property
    def chemicalsubstance_history(self):
        return self._chemicalsubstance_history if self._chemicalsubstance_history != [] else [self.chemicalsubstance]
    @property
    def _currentTg(self):
        """returns the current Tg if its defined, if not None"""
        return check_units(getattr(self, "Tg", None))[0]
    @property
    def _currentporosity(self):
        """returns the current porosity if its defined, if not 0"""
        return getattr(self, "porosity", 0)
    @property
    def _currentcrystallinity(self):
        """returns the crystallinity at the current temperature (if defined)"""
        return getattr(self, "crystallinity", lambda T=None: 0)(T=None)
    @property
    def Tg_history(self):
        return self._Tg_history if self._Tg_history != [] else [self._currentTg]
    @property
    def porosity_history(self):
        return self._porosity_history if self._porosity_history != [] else [self._currentporosity]
    @property
    def crystallinity_history(self):
        return self._crystallinity_history if self._crystallinity_history != [] else [self._currentcrystallinity]
    @property
    def layerclass(self): return type(self).__name__
    @property
    def name(self): return self._name
    @property
    def type(self): return self._type
    @property
    def material(self): return self._material
    @property
    def code(self): return self._code
    @property
    def l(self): return self._l if not self.hasllink else self.llink.getfull(self._l)
    @property
    def D(self):
        Dtmp = None
        if self.Dmodel == "default": # default behavior
            Dtmp = self._compute_Dmodel()
        elif callable(self.Dmodel): # user override
            Dtmp = self.Dmodel()
        if Dtmp is not None:
            Dtmp = np.full_like(self._D, Dtmp,dtype=np.float64)
            if self.hasDlink:
                return self.Dlink.getfull(Dtmp) # substitution rules are applied as defined in Dlink
            else:
                return Dtmp
        return self._D if not self.hasDlink else self.Dlink.getfull(self._D)
    @property
    def k(self):
        ktmp = None
        if self.kmodel == "default": # default behavior
            ktmp = self._compute_kmodel()
        elif callable(self.kmodel): # user override
            ktmp = self.kmodel()
        if ktmp is not None:
            ktmp = np.full_like(self._k, ktmp,dtype=np.float64)
            if self.hasklink:
                return self.klink.getfull(ktmp) # substitution rules are applied as defined in klink
            else:
                return ktmp
        return self._k if not self.hasklink else self.klink.getfull(self._k)
    @property
    def C0(self): return self._C0 if not self.hasC0link else self.COlink.getfull(self._C0)
    @property
    def rho(self): return self._rho
    @property
    def T(self): return self._T if not self.hasTlink else self.Tlink.getfull(self._T)
    @property
    def TK(self): return self._T+T0K
    @property
    def lunit(self): return self._lunit
    @property
    def Dunit(self): return self._Dunit
    @property
    def kunit(self): return self._kunit
    @property
    def Cunit(self): return self._Cunit
    @property
    def rhounit(self): return self._rhounit
    @property
    def Tunit(self): return self._Tunit
    @property
    def TKunit(self): return "K"
    @property
    def n(self): return self._nlayer
    @property
    def nmesh(self): return self._nmesh
    @property
    def nmeshmin(self): return self._nmeshmin
    @property
    def resistance(self): return self.l*self.k/self.D
    @property
    def permeability(self): return self.D/(self.l*self.k)
    @property
    def lag(self): return self.l**2/(6*self.D)
    @property
    def pressure(self): return self.k*self.C0
    @property
    def thickness(self): return sum(self.l)
    @property
    def concentration(self): return sum(self.l*self.C0)/self.thickness
    @property
    def relative_thickness(self): return self.l/self.thickness
    @property
    def relative_resistance(self): return self.resistance/sum(self.resistance)
    @property
    def rank(self): return (self.n-np.argsort(np.array(self.resistance))).tolist()
    @property
    def referencelayer(self): return np.argmax(self.resistance)
    @property
    def lreferencelayer(self): return self.l[self.referencelayer]
    @property
    def Foscale(self): return self.D[self.referencelayer]/self.lreferencelayer**2

    # substance/solute/migrant/chemical (of class migrant or None)
    @property
    def substance(self): return self._substance
    @property
    def migrant(self): return self.substance # alias/synonym of substance
    @property
    def solute(self): return self.substance # alias/synonym of substance
    @property
    def chemical(self): return self.substance # alias/synonym of substance
    # medium (of class foodlayer or None)
    @property
    def medium(self): return self._medium

    # Dmodel and kmodel returned as properties (they are lambda functions)
    # Note about the implementation: They are attributes that remain None or a callable function
    # polymer and mass are udpdated on the fly (the code loops over all layers)
    @property
    def Dmodel(self):
        return self._Dmodel
    @Dmodel.setter
    def Dmodel(self,value):
        if value is None or callable(value):
            self._Dmodel = value
        else:
            raise ValueError("Dmodel must be None or a callable function")
    @property
    def _compute_Dmodel(self):
        """Return a callable function that evaluates D with updated parameters."""
        if not isinstance(self._substance,migrant) or self._substance.Deval() is None:
            return lambda **kwargs: None  # Return a function that always returns None
        template = self._substance.Dtemplate.copy()
        template.update()
        def func(**kwargs):
            D = np.empty_like(self._D)
            for (i,),T in np.ndenumerate(self.T.ravel()): # loop over all layers via T
                template.update(polymer=self.layerclass_history[i],
                                Tg=self.Tg_history[i],
                                T=T) # updated layer properties
                # test if an alternative mode is applicable
                alt_Dclass = self._substance.suggest_alt_Dclass(self,index=i,RaiseError=False,RaiseWarning=False,**template)
                if alt_Dclass is None:
                    # run the default model with inherited eventual user parameters
                    result = self._substance.D.evaluate(**dict(template, **kwargs))
                else:
                    # run the alternative model with inherited eventual user parameters
                    result = alt_Dclass.evaluate(**dict(template, **kwargs))
                if isinstance(result, np.ndarray) and result.size == 1:
                    D[i] = result.item()  # Extract scalar safely
                else:
                    D[i] = result  # Assume it's already a scalar
            return D
        return func # we return a callable function not a value

    # polarity index and molar volume are updated on the fly
    @property
    def kmodel(self):
        return self._kmodel
    @kmodel.setter
    def kmodel(self,value):
        if value is None or callable(value):
            self._kmodel = value
        else:
            raise ValueError("kmodel must be None or a callable function")
    @property
    def _compute_kmodel(self):
        """Return a callable function that evaluates k with updated parameters."""
        if not isinstance(self._substance,migrant) or self._substance.keval() is None:
            return lambda **kwargs: None  # Return a function that always returns None
        template = self._substance.ktemplate.copy()
        # add solute (i) properties: Pi and Vi have been set by loadpubchem already
        template.update(ispolymer = True)
        def func(**kwargs):
            k = np.full_like(self._k,self._k,dtype=np.float64)
            for (i,),T in np.ndenumerate(self.T.ravel()): # loop over all layers via T
                if not self.ispolymer_history[i]: # k can be evaluated only in polymes via FH theory
                    continue # we keep the existing k value
                # add/update monomer properties + porosity and crystallinity of the polymer
                monomer = migrant(self.chemicalsubstance_history[i])
                template.update(Pk = monomer.polarityindex,
                                Vk = monomer.molarvolumeMiller,
                                crystallinity = self.crystallinity_history[i],
                                porosity = self.porosity_history[i])
                # inherit eventual user parameters
                result = self._substance.k.evaluate(**dict(template, **kwargs))
                if isinstance(result, np.ndarray) and result.size == 1:
                    k[i] = result.item()  # Extract scalar safely
                else:
                    k[i] = result  # Assume it's already a scalar
            return k
        return func # we return a callable function not a value


    @property
    def hasDmodel(self):
        """Returns True if a Dmodel has been defined"""
        if hasattr(self, "_compute_Dmodel"):
            if self._compute_Dmodel() is not None:
                return True
            elif callable(self.Dmodel):
                return self.Dmodel() is not None
        return False

    def _currentDmodel(self,index=0):
        """Returns the name of Dmodel used in the ith layer"""
        if self.hasDmodel:
            default = self._substance.D.__name__
            altDmodel = self._substance.suggest_alt_Dmodel(self,index,RaiseError=False)
            if altDmodel is None:
                return default
            return altDmodel if self._substance.check_alt_propclass(altDmodel) else altDmodel+"-->"+default
        return ""

    @property
    def haskmodel(self):
        """Returns True if a kmodel has been defined"""
        if hasattr(self, "_compute_kmodel"):
            if self._compute_kmodel() is not None:
                return True
            elif callable(self.kmodel):
                return self.kmodel() is not None
        return False


    # --------------------------------------------------------------------
    # comparators based resistance
    # --------------------------------------------------------------------
    def __eq__(self, o):
        value1 = self.resistance if self._nlayer>1 else self.resistance[0]
        if isinstance(o,layer):
            value2 = o.resistance if o._nlayer>1 else o.resistance[0]
        else:
            value2 = o
        return value1==value2

    def __ne__(self, o):
        value1 = self.resistance if self._nlayer>1 else self.resistance[0]
        if isinstance(o,layer):
            value2 = o.resistance if o._nlayer>1 else o.resistance[0]
        else:
            value2 = o
        return value1!=value2

    def __lt__(self, o):
        value1 = self.resistance if self._nlayer>1 else self.resistance[0]
        if isinstance(o,layer):
            value2 = o.resistance if o._nlayer>1 else o.resistance[0]
        else:
            value2 = o
        return value1<value2

    def __gt__(self, o):
        value1 = self.resistance if self._nlayer>1 else self.resistance[0]
        if isinstance(o,layer):
            value2 = o.resistance if o._nlayer>1 else o.resistance[0]
        else:
            value2 = o
        return value1>value2

    def __le__(self, o):
        value1 = self.resistance if self._nlayer>1 else self.resistance[0]
        if isinstance(o,layer):
            value2 = o.resistance if o._nlayer>1 else o.resistance[0]
        else:
            value2 = o
        return value1<=value2

    def __ge__(self, o):
        value1 = self.resistance if self._nlayer>1 else self.resistance[0]
        if isinstance(o,layer):
            value2 = o.resistance if o._nlayer>1 else o.resistance[0]
        else:
            value2 = o
        return value1>=value2


    # --------------------------------------------------------------------
    # Generates mesh
    # --------------------------------------------------------------------
    def mesh(self,nmesh=None,nmeshmin=None):
        """ nmesh() generates mesh based on nmesh and nmeshmin, nmesh(nmesh=value,nmeshmin=value) """
        if nmesh==None: nmesh = self.nmesh
        if nmeshmin==None: nmeshmin = self.nmeshmin
        if nmeshmin>nmesh: nmeshmin,nmesh = nmesh, nmeshmin
        # X = mesh distribution (number of nodes per layer)
        X = np.ones(self._nlayer)
        for i in range(1,self._nlayer):
           X[i] = X[i-1]*(self.permeability[i-1]*self.l[i])/(self.permeability[i]*self.l[i-1])
        X = np.maximum(nmeshmin,np.ceil(nmesh*X/sum(X)))
        X = np.round((X/sum(X))*nmesh).astype(int)
        # do the mesh
        x0 = 0
        mymesh = []
        for i in range(self._nlayer):
            mymesh.append(mesh(self.l[i]/self.l[self.referencelayer],X[i],x0=x0,index=i))
            x0 += self.l[i]
        return mymesh

    # --------------------------------------------------------------------
    # Setter methods and tools to validate inputs checknumvalue and checktextvalue
    # --------------------------------------------------------------------
    @physicalstate.setter
    def physicalstate(self,value):
        if value not in ("solid","liquid","gas","supercritical"):
            raise ValueError(f"physicalstate must be solid/liduid/gas/supercritical and not {value}")
        self._physicalstate = value
    @chemicalclass.setter
    def chemicalclass(self,value):
        if value not in ("polymer","other"):
            raise ValueError(f"chemicalclass must be polymer/oher and not {value}")
        self._chemicalclass= value
    @chemicalsubstance.setter
    def chemicalsubstance(self,value):
        if not isinstance(value,str):
            raise ValueError("chemicalsubtance must be str not a {type(value).__name__}")
        self._chemicalsubstance= value
    @polarityindex.setter
    def polarityindex(self,value):
        if not isinstance(value,(float,int)):
            raise ValueError("polarity index must be float not a {type(value).__name__}")
        self._polarityindex= value

    def checknumvalue(self,value,ExpectedUnits=None):
        """ returns a validate value to set properties """
        if isinstance(value,tuple):
            value = check_units(value,ExpectedUnits=ExpectedUnits)[0]
        if isinstance(value,int): value = float(value)
        if isinstance(value,float): value = np.array([value])
        if isinstance(value,list): value = np.array(value)
        if isinstance(value,np.ndarray) and np.ndim(value)==0:
            value = np.atleast_1d(value)
        if len(value)>self._nlayer:
            value = value[:self._nlayer]
            if self.verbosity>1 and self.verbose:
                print('dimension mismatch, the extra value(s) has been removed')
        elif len(value)<self._nlayer:
            value = np.concatenate((value,value[-1:]*np.ones(self._nlayer-len(value))))
            if self.verbosity>1 and self.verbose:
                print('dimension mismatch, the last value has been repeated')
        return value

    def checktextvalue(self,value):
        """ returns a validate value to set properties """
        if not isinstance(value,list): value = [value]
        if len(value)>self._nlayer:
            value = value[:self._nlayer]
            if self.verbosity>1 and self.verbose:
                print('dimension mismatch, the extra entry(ies) has been removed')
        elif len(value)<self._nlayer:
            value = value + value[-1:]*(self._nlayer-len(value))
            if self.verbosity>1 and self.verbose:
                print('dimension mismatch, the last entry has been repeated')
        return value

    @l.setter
    def l(self,value): self._l =self.checknumvalue(value,layer._defaults["lunit"])
    @D.setter
    def D(self,value): self._D=self.checknumvalue(value,layer._defaults["Dunit"])
    @k.setter
    def k(self,value): self._k =self.checknumvalue(value,layer._defaults["kunit"])
    @C0.setter
    def C0(self,value): self._C0 =self.checknumvalue(value,layer._defaults["Cunit"])
    @rho.setter
    def rho(self,value): self._rho =self.checknumvalue(value,layer._defaults["rhounit"])
    @T.setter
    def T(self,value): self._T =self.checknumvalue(value,layer._defaults["Tunit"])
    @name.setter
    def name(self,value): self._name =self.checktextvalue(value)
    @type.setter
    def type(self,value): self._type =self.checktextvalue(value)
    @material.setter
    def material(self,value): self._material =self.checktextvalue(value)
    @nmesh.setter
    def nmesh(self,value): self._nmesh = max(value,self._nlayer*self._nmeshmin)
    @nmeshmin.setter
    def nmeshmin(self,value): self._nmeshmin = max(value,round(self._nmesh/(2*self._nlayer)))
    @substance.setter
    def substance(self,value):
        if isinstance(value,str):
            value = migrant(value)
        if not isinstance(value,migrant) and value is not None:
            raise TypeError(f"value must be a migrant not a {type(value).__name__}")
        self._substance = value
    @migrant.setter
    def migrant(self,value):
        self.substance = value
    @chemical.setter
    def chemical(self,value):
        self.substance = value
    @solute.setter
    def solute(self,value):
        self.substance = value
    @medium.setter
    def medium(self,value):
        from patankar.food import foodlayer
        if not isinstance(value,foodlayer):
            raise TypeError(f"value must be a foodlayer not a {type(value).__name__}")
        self._medium = value

    # --------------------------------------------------------------------
    #  getter and setter for links: Dlink, klink, C0link, Tlink, llink
    # --------------------------------------------------------------------
    @property
    def Dlink(self):
        """Getter for Dlink"""
        return self._Dlink
    @Dlink.setter
    def Dlink(self, value):
        """Setter for Dlink"""
        self._Dlink = self._initialize_link(value, "D")
        if isinstance(value,layerLink): value._maxlength = self.n
    @property
    def klink(self):
        """Getter for klink"""
        return self._klink
    @klink.setter
    def klink(self, value):
        """Setter for klink"""
        self._klink = self._initialize_link(value, "k")
        if isinstance(value,layerLink): value._maxlength = self.n
    @property
    def C0link(self):
        """Getter for C0link"""
        return self._C0link
    @C0link.setter
    def C0link(self, value):
        """Setter for C0link"""
        self._C0link = self._initialize_link(value, "C0")
        if isinstance(value,layerLink): value._maxlength = self.n
    @property
    def Tlink(self):
        """Getter for Tlink"""
        return self._Tlink
    @Tlink.setter
    def Tlink(self, value):
        """Setter for Tlink"""
        self._Tlink = self._initialize_link(value, "T")
        if isinstance(value,layerLink): value._maxlength = self.n
    @property
    def llink(self):
        """Getter for llink"""
        return self._llink
    @llink.setter
    def llink(self, value):
        """Setter for llink"""
        self._llink = self._initialize_link(value, "l")
        if isinstance(value,layerLink): value._maxlength = self.n
    @property
    def hasDlink(self):
        """Returns True if Dlink is defined"""
        return self.Dlink is not None
    @property
    def hasklink(self):
        """Returns True if klink is defined"""
        return self.klink is not None
    @property
    def hasC0link(self):
        """Returns True if C0link is defined"""
        return self.C0link is not None
    @property
    def hasTlink(self):
        """Returns True if Tlink is defined"""
        return self.Tlink is not None
    @property
    def hasllink(self):
        """Returns True if llink is defined"""
        return self.llink is not None

    # --------------------------------------------------------------------
    # returned LaTeX-formated properties
    # --------------------------------------------------------------------
    def Dlatex(self, numdigits=4, units=r"\mathrm{m^2 \cdot s^{-1}}",prefix="D=",mathmode="$"):
        """Returns diffusivity values (D) formatted in LaTeX scientific notation."""
        return [format_scientific_latex(D, numdigits, units, prefix,mathmode) for D in self.D]

    def klatex(self, numdigits=4, units="a.u.",prefix="k=",mathmode="$"):
        """Returns Henry-like values (k) formatted in LaTeX scientific notation."""
        return [format_scientific_latex(k, numdigits, units, prefix,mathmode) for k in self.k]

    def llatex(self, numdigits=4, units="m",prefix="l=",mathmode="$"):
        """Returns thickness values (k) formatted in LaTeX scientific notation."""
        return [format_scientific_latex(l, numdigits, units, prefix,mathmode) for l in self.l]

    def C0latex(self, numdigits=4, units="a.u.",prefix="C0=",mathmode="$"):
        """Returns Initial Concentratoin values (C0) formatted in LaTeX scientific notation."""
        return [format_scientific_latex(c, numdigits, units, prefix,mathmode) for c in self.C0]

    # --------------------------------------------------------------------
    # hash methods (assembly and layer-by-layer)
    # note that list needs to be converted into tuples to be hashed
    # --------------------------------------------------------------------
    def __hash__(self):
        """ hash layer-object (assembly) method """
        return hash((tuple(self._name),
                     tuple(self._type),
                     tuple(self._material),
                     tuple(self._l),
                     tuple(self._D),
                     tuple(self.k),
                     tuple(self._C0),
                     tuple(self._rho)))

    # layer-by-layer @property = decoration to consider it
    # as a property instead of a method/attribute
    # comprehension for n in range(self._nlayer) applies it to all layers
    @property
    def hashlayer(self):
        """ hash layer (layer-by-layer) method """
        return [hash((self._name[n],
                      self._type[n],
                      self._material[n],
                      self._l[n],
                      self._D[n],
                      self.k[n],
                      self._C0[n],
                      self._rho[n]))
                for n in range(self._nlayer)
                ]


    # --------------------------------------------------------------------
    # repr method (since the getter are defined, the '_' is dropped)
    # --------------------------------------------------------------------
    # density and temperature are not shown
    def __repr__(self):
        """ disp method """
        print("\n[%s version=%0.4g, contact=%s]" % (self.__description,self.__version,self.__contact))
        if self._nlayer==0:
            print("empty %s" % (self.__description))
        else:
            hasDmodel, haskmodel = self.hasDmodel, self.haskmodel
            hasDlink, hasklink, hasC0link, hasTlink, hasllink = self.hasDlink, self.hasklink, self.hasC0link, self.hasTlink, self.hasllink
            properties_hasmodel = {"l":False,"D":hasDmodel,"k":haskmodel,"C0":False}
            properties_haslink = {"l":hasllink,"D":hasDlink,"k":hasklink,"C0":hasC0link,"T":hasTlink}
            if hasDmodel or haskmodel:
                properties_hasmodel["T"] = False
            fmtval = '%10s: '+self._printformat+" [%s]"
            fmtstr = '%10s= %s'
            if self._nlayer==1:
                print(f'monolayer of {self.__description}:')
            else:
                print(f'{self._nlayer}-multilayer of {self.__description}:')
            for n in range(1,self._nlayer+1):
                modelinfo = {
                    "D": f"{self._currentDmodel(n-1)}({self.layerclass_history[n-1]},{self._substance},T={float(self.T[0])} {self.Tunit})" if hasDmodel else "",
                    "k": f"{self._substance.k.__name__}(<{self.chemicalsubstance_history[n-1]}>,{self._substance})" if haskmodel else "",
                    }
                print('-- [ layer %d of %d ] ---------- barrier rank=%d --------------'
                      % (n,self._nlayer,self.rank[n-1]))
                # generic properties
                for p in ["name","type","material","code"]:
                    v = getattr(self,p)
                    print('%10s: "%s"' % (p,v[n-1]),flush=True)
                # polymer properties if relevant
                if haskmodel:
                    print('%10s: %s' % ("crystal",self.crystallinity_history[n-1]),flush=True)
                if hasDmodel and self.Tg_history[n-1] is not None:
                    Tg = self.Tg_history[n-1] # Tg,Tgu = check_units(self.Tg_history[n-1])
                    Tg = Tg.item() if isinstance(Tg, np.ndarray) else Tg
                    print('%10s: %s' % ("Tg",f"{Tg} [{self.Tunit}]"),flush=True)
                # numeric properties
                for p in properties_hasmodel.keys():
                    v = getattr(self,p)                 # value
                    vunit = getattr(self,p[0]+"unit")   # value unit
                    print(fmtval % (p,v[n-1],vunit),flush=True)
                    isoverridenbylink = False
                    if properties_haslink[p]:
                        isoverridenbylink = not np.isnan(getattr(self,p+"link").get(n-1))
                    if isoverridenbylink:
                        print(fmtstr % ("",f"value controlled by {p}link[{n-1}] (external)"),flush=True)
                    elif properties_hasmodel[p]:
                        print(fmtstr % ("",modelinfo[p]),flush=True)
        return str(self)

    def __str__(self):
        """Formatted string representation of layer"""
        all_identical = len(set(self.layerclass_history)) == 1
        cls = self.__class__.__name__ if all_identical else "multilayer"
        return f"<{cls} with {self.n} layer{'s' if self.n>1 else ''}: {self.name}>"

    # --------------------------------------------------------------------
    # Returns the equivalent dictionary from an object for debugging
    # --------------------------------------------------------------------
    def _todict(self):
        """ returns the equivalent dictionary from an object """
        return dict((key, getattr(self, key)) for key in dir(self) if key not in dir(self.__class__))
    # --------------------------------------------------------------------

    # --------------------------------------------------------------------
    # Simplify layers by collecting similar ones
    # --------------------------------------------------------------------
    def simplify(self):
        """ merge continuous layers of the same type """
        nlayer = self._nlayer
        if nlayer>1:
           res = self[0]
           ires = 0
           ireshash = res.hashlayer[0]
           for i in range(1,nlayer):
               if self.hashlayer[i]==ireshash:
                   res.l[ires] = res.l[ires]+self.l[i]
               else:
                   res = res + self[i]
                   ires = ires+1
                   ireshash = self.hashlayer[i]
        else:
             res = self.copy()
        return res

    # --------------------------------------------------------------------
    # Split layers into a tuple
    # --------------------------------------------------------------------
    def split(self):
        """ split layers """
        out = ()
        if self._nlayer>0:
            for i in range(self._nlayer):
                out = out + (self[i],) # (,) special syntax for tuple singleton
        return out

    # --------------------------------------------------------------------
    # deepcopy
    # --------------------------------------------------------------------
    def copy(self,**kwargs):
        """
        Creates a deep copy of the current layer instance.

        Returns:
        - layer: A new layer instance identical to the original.
        """
        return duplicate(self).update(**kwargs)

    # --------------------------------------------------------------------
    # update contact conditions from a foodphysics instance (or do the reverse)
    # material << medium
    # material@medium
    # --------------------------------------------------------------------
    def _from(self,medium=None):
        """Propagates contact conditions from food instance"""
        from patankar.food import foodphysics, foodlayer
        if not isinstance(medium,foodphysics):
            raise TypeError(f"medium must be a foodphysics, foodlayer not a {type(medium).__name__}")
        if not hasattr(medium, "contacttemperature"):
            medium.contacttemperature = self.T[0]
        T = medium.get_param("contacttemperature",40,acceptNone=False)
        self.T = np.full_like(self.T,T,dtype=np.float64)
        if medium.substance is not None:
            self.substance = medium.substance
        else:
            medium.substance = self.substance # do the reverse if substance is not defined in medium
        # inherit fully medium only if it is a foodlayer (foodphysics is too restrictive)
        if isinstance(medium,foodlayer):
            self.medium = medium

    # overload operator <<
    def __lshift__(self, medium):
        """Overloads << to propagate contact conditions from food."""
        self._from(medium)
    # overload operator @ (same as <<)
    def __matmul__(self, medium):
        """Overloads @ to propagate contact conditions from food."""
        self._from(medium)


    # --------------------------------------------------------------------
    # Inheritance registration mechanism associated with food >> layer
    # It is used by food, not by layer (please refer to food.py).
    # Note that layer >> food means mass transfer simulation
    # --------------------------------------------------------------------
    def acknowledge(self, what=None, category=None):
        """
        Register inherited properties under a given category.

        Parameters:
        -----------
        what : str or list of str or a set
            The properties or attributes that have been inherited.
        category : str
            The category under which the properties are grouped.
        """
        if category is None or what is None:
            raise ValueError("Both 'what' and 'category' must be provided.")
        if isinstance(what, str):
            what = {what}  # Convert string to a set
        elif isinstance(what, list):
            what = set(what)  # Convert list to a set for uniqueness
        elif not isinstance(what,set):
            raise TypeError("'what' must be a string, a list, or a set of strings.")
        if category not in self._hasbeeninherited:
            self._hasbeeninherited[category] = set()
        self._hasbeeninherited[category].update(what)

    # --------------------------------------------------------------------
    # migration simulation overloaded as sim = layer >> food
    # using layer >> food without output works also.
    # The result is stored in food.lastsimulation
    # --------------------------------------------------------------------
    def contact(self,medium,**kwargs):
        """alias to migration method"""
        return self.migration(medium,**kwargs)

    def migration(self,medium=None,**kwargs):
        """interface to simulation engine: senspantankar"""
        from patankar.food import foodphysics
        from patankar.migration import senspatankar
        if medium is None:
            medium = self.medium
        if not isinstance(medium,foodphysics):
            raise TypeError(f"medium must be a foodphysics not a {type(medium).__name__}")
        sim = senspatankar(self,medium,**kwargs)
        medium.lastsimulation = sim # store the last simulation result in medium
        medium.lastinput = self # store the last input (self)
        sim.savestate(self,medium) # store store the inputs in sim for chaining
        return sim

    # overloading operation >>
    def __rshift__(self, medium):
        """
            Overloads >> to propagate migration to food.
                layer >> food --> simulation --> food (use food.lastsimulation to see it)
                layer >> condition --> simulation --> condition
        """
        from patankar.food import foodphysics
        if not isinstance(medium,foodphysics):
            raise TypeError(f"medium must be a foodphysics object not a {type(medium).__name__}")
        return self.contact(medium)

    # overloading in
    def __rmod__(self,other):
        """
            Overload in to enable: susbtance % layer --> layer
        """
        from patankar.loadpubchem import migrant
        if not isinstance(other,migrant):
            raise TypeError(f"other must a migrant and not a {type(other).__name__}")
        self.update(substance=other)
        return self

    # --------------------------------------------------------------------
    # Safe update method
    # --------------------------------------------------------------------
    def update(self, **kwargs):
        """
        Update layer parameters following strict validation rules.

        Rules:
        1) key should be listed in self._defaults
        2) for some keys, synonyms are acceptable as reported in self._synonyms
        3) values cannot be None if they were not None in _defaults
        4) values should be str if they were initially str, idem with bool
        5) values which were numeric (int, float, np.ndarray) should remain numeric.
        6) lists are acceptable as numeric arrays
        7) all numerical (float, np.ndarray, list) except int must be converted into numpy arrays.
           Values which were int in _defaults must remain int and an error should be raised
           if a float value is proposed.
        8) keys listed in _parametersWithUnits can be assigned with tuples (value, "unit").
           They will be converted automatically with check_units(value).
        9) for parameters with a default value None, any value is acceptable
        10) A clear error message should be displayed for any bad value showing the
            current value of the parameter and its default value.
        """

        if not kwargs:  # shortcut
            return self # for chaining

        param_counts = {key: 0 for key in self._defaults}  # Track how many times each param is set

        def resolve_key(key):
            """Resolve key considering synonyms and check for duplicates."""
            for main_key, synonyms in self._synonyms.items():
                if key == main_key or key in synonyms:
                    param_counts[main_key] += 1
                    return main_key
            param_counts[key] += 1
            return key

        def validate_value(key, value):
            """Validate and process the value according to the rules."""
            default_value = self._defaults[key]

            # Rule 3: values cannot be None if they were not None in _defaults
            if value is None and default_value is not None:
                raise ValueError(f"Invalid value for '{key}': None is not allowed. "
                                 f"Current: {getattr(self, key)}, Default: {default_value}")

            # Rule 9: If default is None, any value is acceptable
            if default_value is None:
                return value

            # Rule 4 & 5: Ensure type consistency (str, bool, or numeric types)
            if isinstance(default_value, str) and not isinstance(value, str):
                raise TypeError(f"Invalid type for '{key}': Expected str, got {type(value).__name__}. "
                                f"Current: {getattr(self, key)}, Default: {default_value}")
            if isinstance(default_value, bool) and not isinstance(value, bool):
                raise TypeError(f"Invalid type for '{key}': Expected bool, got {type(value).__name__}. "
                                f"Current: {getattr(self, key)}, Default: {default_value}")

            # Rule 6 & 7: Convert numeric types properly
            if isinstance(default_value, (int, float, np.ndarray)):
                if isinstance(value, list):
                    value = np.array(value)

                if isinstance(default_value, int):
                    if isinstance(value, float) or (isinstance(value, np.ndarray) and np.issubdtype(value.dtype, np.floating)):
                        raise TypeError(f"Invalid type for '{key}': Expected integer, got float. "
                                        f"Current: {getattr(self, key)}, Default: {default_value}")
                    if isinstance(value, (int, np.integer)):
                        return int(value)  # Ensure it remains an int
                    raise TypeError(f"Invalid type for '{key}': Expected integer, got {type(value).__name__}. "
                                    f"Current: {getattr(self, key)}, Default: {default_value}")

                if isinstance(value, (int, float, list, np.ndarray)):
                    return np.array(value, dtype=float)  # Convert everything to np.array for floats

                raise TypeError(f"Invalid type for '{key}': Expected numeric, got {type(value).__name__}. "
                                f"Current: {getattr(self, key)}, Default: {default_value}")

            # Rule 8: Convert units if applicable
            if key in self._parametersWithUnits and isinstance(value, tuple):
                value, unit = value
                converted_value, _ = check_units((value, unit), ExpectedUnits=self._parametersWithUnits[key])
                return converted_value

            return value

        # Apply updates while tracking parameter occurrences
        for key, value in kwargs.items():
            resolved_key = resolve_key(key)

            if resolved_key not in self._defaults:
                raise KeyError(f"Invalid key '{key}'. Allowed keys: {list(self._defaults.keys())}.")

            try:
                validated_value = validate_value(resolved_key, value)
                setattr(self, resolved_key, validated_value)
            except (TypeError, ValueError) as e:
                raise ValueError(f"Error updating '{key}': {e}")

        # Ensure that no parameter was set multiple times due to synonyms
        duplicate_keys = [k for k, v in param_counts.items() if v > 1]
        if duplicate_keys:
            raise ValueError(f"Duplicate assignment detected for parameters: {duplicate_keys}. "
                             "Use only one synonym per parameter.")

        return self # to enable chaining

    # Basic tool for debugging
    # --------------------------------------------------------------------
    # STRUCT method - returns the equivalent dictionary from an object
    # --------------------------------------------------------------------
    def struct(self):
        """ returns the equivalent dictionary from an object """
        return dict((key, getattr(self, key)) for key in dir(self) if key not in dir(self.__class__))

# %% Mesh class
# Mesh class
# =======================
class mesh():
    """ simple nodes class for finite-volume methods """
    def __init__(self,l,n,x0=0,index=None):
       self.x0 = x0
       self.l = l
       self.n = n
       de = dw = l/(2*n)
       self.de = np.ones(n)*de
       self.dw = np.ones(n)*dw
       self.xmesh = np.linspace(0+dw,l-de,n) # nodes positions
       self.w = self.xmesh - dw
       self.e = self.xmesh + de
       self.index = np.full(n, int(index), dtype=np.int32)

    def __repr__(self):
        print(f"-- mesh object (layer index={self.index[0]}) --")
        print("%25s = %0.4g" % ("start at x0", self.x0))
        print("%25s = %0.4g" % ("domain length l", self.l))
        print("%25s = %0.4g" % ("number of nodes n", self.n))
        print("%25s = %0.4g" % ("dw", self.dw[0]))
        print("%25s = %0.4g" % ("de", self.de[0]))
        return "mesh%d=[%0.4g %0.4g]" % \
            (self.n,self.x0+self.xmesh[0],self.x0+self.xmesh[-1])


# %% Material classes
"""
=======================================================
Child classes derived from layer
this section can be extended to define specific layers
    * polymer
    * ink
    * air
    * paper and board

These classes are more flexible than the parent class layer.
They can include temperature dependence, refined tunning, etc.

Properties taken from
    * linear thermal expansoin
    https://omnexus.specialchem.com/polymer-properties/properties/coefficient-of-linear-thermal-expansion

Once the layers are incorporated in a multilayer structure,
they loose their original subclass and become only an object
layer. These subclasses are therefore useful to refine the
properties of each layer before standarizing them.

Polarity index is used as an helper to set Henri-like coefficients.
A common scale for polarity index for solvents is from 0 to 10:
- 0-3: Non-polar solvents (e.g., hexane)
- 4-6: Moderately polar solvents (e.g., acetone)
- 7-10: Polar solvents (e.g., water)
We consider that polymers are solid solvents.
=========================================================
"""

# <<<<<<<<<<<<<<<<<<<<<<< P O L Y O L E F I N S >>>>>>>>>>>>>>>>>>>>>>

# <-- LDPE polymer ---------------------------------->
class LDPE(layer):
    """  extended pantankar.layer for low-density polyethylene LDPE  """
    _chemicalsubstance = "ethylene" # monomer for polymers
    _polarityindex = 1.0  # Very non-polar (typical for polyolefins)
    def __init__(self,l=100e-6,D=1e-12,T=None,
                 k=None,C0=None,lunit=None,Dunit=None,kunit=None,Cunit=None,
                 layername="layer in LDPE",**extra):
        """ LDPE layer constructor """
        super().__init__(
                       l=l,D=D,k=k,C0=C0, T=T,
                       lunit=lunit,Dunit=Dunit,kunit=kunit,Cunit=Cunit,
                       layername=layername,
                       layertype="polymer", # set by default at inititialization
                       layermaterial="low-density polyethylene",
                       layercode="LDPE",
                       **extra
                       )
    def density(self,T=None):
        """ density of LDPE: density(T in K) """
        T = self.T[0] if T is None else check_units(T,None,"degC")[0]
        return 920 *(1-3*(T-layer._defaults["Td"])*20e-5),"kg/m**3" # lowest temperature
    @property
    def Tg(self):
        """ glass transition temperature of LDPE """
        return -130,"degC" # lowest temperature
    @property
    def Tm(self):
        """ typical melting temperature of LDPE """
        return 110, "degC"

    def crystallinity(self, T=None):
        """polymer crystallinity"""
        T = self.T[0] if T is None else T
        if self.Tm is None or T > self.Tm[0]:
            return 0
        else:
            return 0.4
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0

# <-- HDPE polymer ---------------------------------->
class HDPE(layer):
    """  extended pantankar.layer for high-density polyethylene HDPE  """
    _chemicalsubstance = "ethylene" # monomer for polymers
    _polarityindex = 2.0 # Non-polar, slightly higher density, similar overall polarity to LDPE
    def __init__(self,l=500e-6,D=1e-13, T=None,
                 k=None,C0=None,lunit=None,Dunit=None,kunit=None,Cunit=None,
                 layername="layer in HDPE",**extra):
        """ HDPE layer constructor """
        layer.__init__(self,
                       l=l,D=D,k=k,C0=C0, T=T,
                       lunit=lunit,Dunit=Dunit,kunit=kunit,Cunit=Cunit,
                       layername=layername,
                       layertype="polymer", # set by default at inititialization
                       layermaterial="high-density polyethylene",
                       layercode="HDPE",
                       **extra
                       )
    def density(self,T=None):
        """ density of HDPE: density(T in K) """
        T = self.T[0] if T is None else check_units(T,None,"degC")[0]
        return 940 *(1-3*(T-layer._defaults["Td"])*11e-5),"kg/m**3" # lowest temperature

    @property
    def Tg(self):
        """ glass transition temperature of HDPE """
        return -100,"degC" # highest temperature

    @property
    def Tm(self):
        """ typical melting temperature of HDPE """
        return 133.5, "degC"
    def crystallinity(self, T=None):
        """polymer crystallinity"""
        T = self.T[0] if T is None else T
        if self.Tm is None or T > self.Tm[0]:
            return 0
        else:
            return 0.8
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0

# <-- LLDPE polymer ---------------------------------->
class LLDPE(layer):
    """ extended pantankar.layer for linear low-density polyethylene LLDPE """
    _chemicalsubstance = "ethylene" # monomer for polymers
    _polarityindex = 1.5 # Similar to LDPE, can be slightly more polar if co-monomer is present
    def __init__(self, l=80e-6, D=1e-12, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in LLDPE",**extra):
        """
        LLDPE layer constructor
        Defaults are set to typical values found in the literature or between
        LDPE/HDPE ones. Adjust them as necessary for your models.
        """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="polymer",
            layermaterial="linear low-density polyethylene",
            layercode="LLDPE",
            **extra
        )
    def density(self, T=None):
        """
        density of LLDPE: density(T in K)
        By default, uses an approximate value between LDPE and HDPE.
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 915 * (1 - 3 * (T - layer._defaults["Td"]) * 15e-5), "kg/m**3"
    @property
    def Tg(self):
        """
        glass transition temperature of LLDPE
        Typically close to LDPE, though slightly higher or lower can be found in the literature.
        """
        return -120, "degC"
    @property
    def Tm(self):
        """ typical melting temperature of LLDPE """
        return 110, "degC"
    def crystallinity(self, T=None):
        """polymer crystallinity"""
        T = self.T[0] if T is None else T
        if self.Tm is None or T > self.Tm[0]:
            return 0
        else:
            return 0.45
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0

# <-- PP polymer ---------------------------------->
class PP(layer):
    _chemicalsubstance = "propylene" # monomer for polymers
    _polarityindex = 1.0  # Among the least polar, similar to PE
    """  extended pantankar.layer for isotactic polypropylene PP  """
    def __init__(self,l=300e-6,D=1e-14, T=None,
                 k=None,C0=None,lunit=None,Dunit=None,kunit=None,Cunit=None,
                 layername="layer in PP",**extra):
        """ PP layer constructor """
        layer.__init__(self,
                       l=l,D=D,k=k,C0=C0, T=T,
                       lunit=lunit,Dunit=Dunit,kunit=kunit,Cunit=Cunit,
                       layername=layername,
                       layertype="polymer", # set by default at inititialization
                       layermaterial="isotactic polypropylene",
                       layercode="PP",
                       **extra
                       )
    def density(self,T=None):
        """ density of PP: density(T in K) """
        T = self.T[0] if T is None else check_units(T,None,"degC")[0]
        return 910 *(1-3*(T-layer._defaults["Td"])*7e-5),"kg/m**3" # lowest temperature
    @property
    def Tg(self):
        """ glass transition temperature of PP """
        return 0,"degC" # highest temperature
    @property
    def Tm(self):
        """ typical melting temperature of PP """
        return 165, "degC"
    def crystallinity(self, T=None):
        """polymer crystallinity"""
        T = self.T[0] if T is None else T
        if self.Tm is None or T > self.Tm[0]:
            return 0
        else:
            return 0.5
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0

# -- PPrubber (atactic polypropylene) ---------------------------------
class PPrubber(layer):
    _chemicalsubstance = "propylene" # monomer for polymers
    _polarityindex = 1.0  # Also very non-polar
    """ extended pantankar.layer for atactic (rubbery) polypropylene PP """
    def __init__(self, l=100e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in PPrubber",**extra):
        """ PPrubber layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="polymer",
            layermaterial="atactic polypropylene",
            layercode="aPP",
            **extra
        )
    def density(self, T=None):
        """
        density of atactic (rubbery) PP: density(T in K)
        Approximate initial density ~900 kg/m^3, linear thermal expansion factor
        can be adjusted.
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 900 * (1 - 3*(T - layer._defaults["Td"]) * 17e-5), "kg/m**3"
    @property
    def Tg(self):
        """ glass transition temperature of atactic/rubbery PP """
        return -20, "degC"
    @property
    def Tm(self):
        """ atactic PP has no well-defined melting temperature """
        return None
    def crystallinity(self, T=None):
        """polymer crystallinity"""
        T = self.T[0] if T is None else T
        return 0
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0

# -- oPP (bioriented polypropylene) ------------------------------------
class oPP(layer):
    """ extended pantankar.layer for bioriented polypropylene oPP """
    _chemicalsubstance = "propylene" # monomer for polymers
    _polarityindex = 1.0   # Non-polar, but oriented film might have slight morphological differences
    def __init__(self, l=40e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in oPP",**extra):
        """ oPP layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="polymer",
            layermaterial="bioriented polypropylene",
            layercode="oPP",
            **extra
        )
    def density(self, T=None):
        """
        density of bioriented PP: density(T in K)
        Typically close to isotactic PP around ~910 kg/m^3.
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 910 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"
    @property
    def Tg(self):
        """ glass transition temperature of bioriented PP """
        return 0, "degC"
    @property
    def Tm(self):
        """ typical melting temperature of bioriented PP """
        return 165, "degC"
    def crystallinity(self, T=None):
        """polymer crystallinity"""
        T = self.T[0] if T is None else T
        if self.Tm is None or T > self.Tm[0]:
            return 0
        else:
            return 0.5
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0


# <<<<<<<<<<<<<<<<<<<<<<< P O L Y A C R Y L A T E S >>>>>>>>>>>>>>>>>>>>>>
# -- PMMA (polymethyl acrylate) -----------------------------------------------
class PMMA(layer):
    """ extended pantankar.layer for polystyrene (PS) """
    _chemicalsubstance = "methyl methacrylate" # monomer for polymers
    _polarityindex = 5.5  #  more polar than polyolefins
    def __init__(self, l=100e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in PMMA",**extra):
        """ PMMA layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="polymer",
            layermaterial="PMMA",
            layercode="PMMA",
            **extra
        )
    def density(self, T=None):
        """
        density of PMMA: ~1180 kg/m^3
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 1180 * (1 - 3*(T - layer._defaults["Td"]) * 7.5e-5), "kg/m**3"
    @property
    def Tg(self):
        """ glass transition temperature of PMMA """
        return 105, "degC"
    @property
    def Tm(self):
        """ polystyrene is largely amorphous """
        return None
    def crystallinity(self, T=None):
        """polymer crystallinity"""
        return 0
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0

# <<<<<<<<<<<<<<<<<<<<<<< P O L Y V I N Y L S >>>>>>>>>>>>>>>>>>>>>>

# -- PS (polystyrene) -----------------------------------------------
class PS(layer):
    """ extended pantankar.layer for polystyrene (general purpose PS) """
    _chemicalsubstance = "styrene" # monomer for polymers
    _polarityindex = 3.0  # Slightly more polar than polyolefins, but still considered relatively non-polar
    def __init__(self, l=100e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in PS",layertype="polymer",layermaterial="polystyrene",layercode="PS",
                 **extra):
        """ PS layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype=layertype,
            layermaterial=layermaterial,
            layercode=layercode,
            **extra
        )
    def density(self, T=None):
        """
        density of PS: ~1050 kg/m^3
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 1050 * (1 - 3*(T - layer._defaults["Td"]) * 5e-5), "kg/m**3"
    @property
    def Tg(self):
        """ glass transition temperature of PS """
        return 100, "degC"
    @property
    def Tm(self):
        """ polystyrene is largely amorphous """
        return None
    def crystallinity(self, T=None):
        """polymer crystallinity"""
        return 0
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0

# -- rPS (polystyrene above Tg) --------------------------------------
class rPS(PS):
    """extended pantankar layer for rubber polystyrene (general purpose PS)"""
    def __init__(self, l=100e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in rPS",layermaterial="polystyrene above Tg",
                 **extra):
        """ PS layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layermaterial=layermaterial,
            **extra
        )

# -- HIPS (high-impact polystyrene) -----------------------------------
class HIPS(layer):
    """ extended pantankar.layer for high-impact polystyrene (HIPS) """
    _chemicalsubstance = "styrene" # monomer for polymers
    _polarityindex = 3.0  # Similar or very close to PS in polarity
    def __init__(self, l=100e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in HIPS", layertype="polymer",
                 layermaterial="high-impact polystyrene", layercode="HIPS", **extra):
        """ HIPS layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype=layertype,
            layermaterial=layermaterial,
            layercode=layercode,
            **extra
        )
    def density(self, T=None):
        """
        density of HIPS: ~1040 kg/m^3
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 1040 * (1 - 3*(T - layer._defaults["Td"]) * 5e-5), "kg/m**3"
    @property
    def Tg(self):
        """ glass transition temperature of HIPS """
        return 95, "degC"
    @property
    def Tm(self):
        """ HIPS is also considered amorphous """
        return None
    def crystallinity(self, T=None):
        """polymer crystallinity"""
        return 0
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0

# -- rHIPS (high-impact polystyrene above Tg) --------------------------
class rHIPS(HIPS):
    """extended pantankar layer for rubber high-impact polystyrene"""
    def __init__(self, l=100e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in rHIPS",layermaterial="high-impact polystyrene above Tg",
                 **extra):
        """ PS layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layermaterial=layermaterial,
            **extra
        )


# -- PBS (assuming a styrene-based polymer) ---------------------------
class SBS(layer):
    _chemicalsubstance = "styrene" # Styrene + butadiene
    _polarityindex = 3.5  # Non-polar but somewhat more interactive than pure PE/PP due to styrene units
    """
    extended pantankar.layer for a styrene-based SBS
    Adjust Tg/density as needed for your scenario.
    """
    def __init__(self, l=100e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in PBS",**extra):
        """ DBS layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="polymer",
            layermaterial="styrene-based polymer SBS",
            layercode="SBS",
            **extra
        )
    def density(self, T=None):
        """
        density of 'DBS': approximate, around ~1030 kg/m^3
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 1030 * (1 - 3*(T - layer._defaults["Td"]) * 5e-5), "kg/m**3"
    @property
    def Tg(self):
        """ glass transition temperature of 'DBS' """
        return 90, "degC"
    @property
    def Tm(self):
        """ SBS typically no well-defined melt """
        return None
    def crystallinity(self, T=None):
        """polymer crystallinity"""
        return 0
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0

# -- rigidPVC ---------------------------------------------------------
class rigidPVC(layer):
    """ extended pantankar.layer for rigid PVC """
    _chemicalsubstance = "vinyl chloride" # monomer for polymers
    _polarityindex = 4.0  # Chlorine substituents give moderate polarity.
    def __init__(self, l=200e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in rigid PVC",**extra):
        """ rigid PVC layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="polymer",
            layermaterial="rigid PVC",
            layercode="PVC",
            **extra
        )
    def density(self, T=None):
        """
        density of rigid PVC: ~1400 kg/m^3
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 1400 * (1 - 3*(T - layer._defaults["Td"]) * 15e-5), "kg/m**3"

    @property
    def Tg(self):
        """ glass transition temperature of rigid PVC """
        return 80, "degC"
    @property
    def Tm(self):
        """ rigid PVC is mostly amorphous """
        return None
    def crystallinity(self, T=None):
        """polymer crystallinity"""
        return 0
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0

# -- plasticizedPVC ---------------------------------------------------
class plasticizedPVC(layer):
    """ extended pantankar.layer for plasticized PVC """
    _chemicalsubstance = "vinyl chloride" # monomer for polymers
    _polarityindex = 4.5  # Plasticizers can slightly change overall polarity/solubility.
    def __init__(self, l=200e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in plasticized PVC",**extra):
        """ plasticized PVC layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="polymer",
            layermaterial="plasticized PVC",
            layercode="pPVC",
            **extra
        )
    def density(self, T=None):
        """
        density of plasticized PVC: ~1300 kg/m^3
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 1300 * (1 - 3*(T - layer._defaults["Td"]) * 15e-5), "kg/m**3"
    @property
    def Tg(self):
        """ glass transition temperature of plasticized PVC """
        return -40, "degC"
    @property
    def Tm(self):
        """ plasticized PVC also amorphous """
        return None
    def crystallinity(self, T=None):
        """polymer crystallinity"""
        return 0
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0

# -- plasticizedPVC ---------------------------------------------------
class PVAc(layer):
    """ extended pantankar.layer for PVAc """
    _chemicalsubstance = "vinyl acetate" # monomer for polymers
    _polarityindex = 7  # polar but it depends on the acetylation rate.
    def __init__(self, l=200e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in PVAc",**extra):
        """ PVAc layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="polymer",
            layermaterial="PVAc",
            layercode="PVAc",
            **extra
        )
    def density(self, T=None):
        """
        density of plasticized PVC: ~1300 kg/m^3
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 1190 * (1 - 3*(T - layer._defaults["Td"]) * 15e-5), "kg/m**3"
    @property
    def Tg(self):
        """ glass transition temperature of plasticized PVC """
        return -40, "degC"
    @property
    def Tm(self):
        """ plasticized PVC also amorphous """
        return None
    def crystallinity(self, T=None):
        """polymer crystallinity"""
        return 0
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0


# <<<<<<<<<<<<<<<<<<<<<<< P O L Y E S T E R S >>>>>>>>>>>>>>>>>>>>>>

# -- gPET (glassy PET, T < 76°C) --------------------------------------
class gPET(layer):
    """ extended pantankar.layer for PET in its glassy state (below ~76°C) """
    _chemicalsubstance = "ethylene terephthalate" # monomer for polymers
    _polarityindex = 5.0  # Polyester with significant dipolar interactions (Ph = phenylene ring).
    def __init__(self, l=200e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in gPET",layertype="polymer",layermaterial="glassy PET", layercode="PET",**extra):
        """ glassy PET layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype=layertype,
            layermaterial=layermaterial,
            layercode=layercode,
            **extra
        )
    def density(self, T=None):
        """
        density of glassy PET: ~1350 kg/m^3
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 1350 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"
    @property
    def Tg(self):
        """ approximate glass transition temperature of PET """
        return 76, "degC"
    @property
    def Tm(self):
        """ approximate melting temperature of PET """
        return 250, "degC"
    def crystallinity(self, T=None):
        """polymer crystallinity"""
        T = self.T[0] if T is None else T
        if self.Tm is None or T > self.Tm[0]:
            return 0
        else:
            return 0.35
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0

## wPET(plasticized/wet PET) -------------------------------------------
class wPET(gPET):
    """ extended pantankar.layer for severaly plasticized PET (Tg ~46°C) """
    def __init__(self, l=200e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in wPET",layermaterial="plasticized PET",**extra):
        """ plasticized PET layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layermaterial=layermaterial,
            **extra
        )

    @property
    def Tg(self):
        """ approximate glass transition temperature of PET """
        return 46, "degC"


# -- rPET (rubbery PET, T > 76°C) --------------------------------------
class rPET(layer):
    """ extended pantankar.layer for PET in its rubbery state (above ~76°C) """
    _chemicalsubstance = "ethylene terephthalate" # monomer for polymers
    _polarityindex = 5.0  # Polyester with significant dipolar interactions (Ph = phenylene ring).
    def __init__(self, l=200e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in rPET",**extra):
        """ rubbery PET layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="polymer",
            layermaterial="rubbery PET",
            layercode="rPET",
            **extra
        )
    def density(self, T=None):
        """
        density of rubbery PET: ~1350 kg/m^3
        but with a different expansion slope possible, if needed
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 1350 * (1 - 3*(T - layer._defaults["Td"]) * 1e-4), "kg/m**3"
    @property
    def Tg(self):
        """ approximate glass transition temperature of PET """
        return 76, "degC"
    @property
    def Tm(self):
        """ approximate melting temperature of PET """
        return 250, "degC"
    def crystallinity(self, T=None):
        """polymer crystallinity"""
        T = self.T[0] if T is None else T
        if self.Tm is None or T > self.Tm[0]:
            return 0
        else:
            return 0.35
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0

# -- PBT --------------------------------------------------------------
class PBT(layer):
    """ extended pantankar.layer for polybutylene terephthalate (PBT) """
    _chemicalsubstance = "Buthylene terephthalate" # monomer for polymers
    _polarityindex = 5.5  # Similar to PET, slight structural differences
    def __init__(self, l=200e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in PBT",**extra):
        """ PBT layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="polymer",
            layermaterial="polybutylene terephthalate",
            layercode="PBT",
            **extra
        )
    def density(self, T=None):
        """
        density of PBT: ~1310 kg/m^3
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 1310 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"
    @property
    def Tg(self):
        """ glass transition temperature of PBT """
        return 40, "degC"
    @property
    def Tm(self):
        """ approximate melting temperature of PBT """
        return 225, "degC"
    def crystallinity(self, T=None):
        """polymer crystallinity"""
        T = self.T[0] if T is None else T
        if self.Tm is None or T > self.Tm[0]:
            return 0
        else:
            return 0.35
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0

# -- PEN --------------------------------------------------------------
class PEN(layer):
    _chemicalsubstance = "Buthylene terephthalate" # monomer for polymers
    _polarityindex = 6  # More aromatic than PET, often better barrier properties
    """ extended pantankar.layer for polyethylene naphthalate (PEN) """
    def __init__(self, l=200e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in PEN",**extra):
        """ PEN layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="polymer",
            layermaterial="polyethylene naphthalate",
            layercode="PEN",
            **extra
        )
    def density(self, T=None):
        """
        density of PEN: ~1330 kg/m^3
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 1330 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"
    @property
    def Tg(self):
        """ glass transition temperature of PEN """
        return 120, "degC"
    @property
    def Tm(self):
        """ approximate melting temperature of PEN """
        return 270, "degC"
    def crystallinity(self, T=None):
        """polymer crystallinity"""
        T = self.T[0] if T is None else T
        if self.Tm is None or T > self.Tm[0]:
            return 0
        else:
            return 0.4
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0


# <<<<<<<<<<<<<<<<<<<<<<< P O L Y A M I D E S >>>>>>>>>>>>>>>>>>>>>>

# -- PA6 --------------------------------------------------------------
class PA6(layer):
    _chemicalsubstance = "caprolactam" # monomer for polymers
    _polarityindex = 7.5  # Strong hydrogen-bonding, thus quite polar.
    """ extended pantankar.layer for polyamide 6 (PA6) """
    def __init__(self, l=200e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in PA6",**extra):
        """ PA6 layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="polymer",
            layermaterial="polyamide 6",
            layercode="PA6",
            **extra
        )
    def density(self, T=None):
        """
        density of PA6: ~1140 kg/m^3
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 1140 * (1 - 3*(T - layer._defaults["Td"]) * 5e-5), "kg/m**3"
    @property
    def Tg(self):
        """ glass transition temperature of PA6 """
        return 50, "degC"
    @property
    def Tm(self):
        """ approximate melting temperature of PA6 """
        return 220, "degC"
    def crystallinity(self, T=None):
        """polymer crystallinity"""
        T = self.T[0] if T is None else T
        if self.Tm is None or T > self.Tm[0]:
            return 0
        else:
            return 0.3
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0


# -- PA66 -------------------------------------------------------------
class PA66(layer):
    _chemicalsubstance = "hexamethylenediamine" # monomer for polymers
    _polarityindex = 7.5  # Similar to PA6, strongly polar with hydrogen bonds.
    """ extended pantankar.layer for polyamide 66 (PA66) """
    def __init__(self, l=200e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in PA66",**extra):
        """ PA66 layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="polymer",
            layermaterial="polyamide 6,6",
            layercode="PA6,6",
            **extra
        )
    def density(self, T=None):
        """
        density of PA66: ~1150 kg/m^3
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 1150 * (1 - 3*(T - layer._defaults["Td"]) * 5e-5), "kg/m**3"
    @property
    def Tg(self):
        """ glass transition temperature of PA66 """
        return 70, "degC"
    @property
    def Tm(self):
        """ approximate melting temperature of PA66 """
        return 255, "degC"

    def crystallinity(self, T=None):
        """polymer crystallinity"""
        T = self.T[0] if T is None else T
        if self.Tm is None or T > self.Tm[0]:
            return 0
        else:
            return 0.35
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0


# <<<<<<<<<<<<<<<<<<<<<<< A D H E S I V E S >>>>>>>>>>>>>>>>>>>>>>

# -- AdhesiveNaturalRubber --------------------------------------------
class AdhesiveNaturalRubber(layer):
    _chemicalsubstance = "cis-1,4-polyisoprene" # monomer for polymers
    _polarityindex = 2  # Mostly non-polar; elasticity from cis-isoprene chains.
    """ extended pantankar.layer for natural rubber adhesives """
    def __init__(self, l=20e-6, D=1e-13, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="adhesive natural rubber",**extra):
        """ constructor for a natural rubber-based adhesive layer """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="adhesive",
            layermaterial="natural rubber adhesive",
            layercode="rubber",
            **extra
        )
    def density(self, T=None):
        """ typical density ~910 kg/m^3, adjust as needed """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 910 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"
    @property
    def Tg(self):
        """ approximate Tg of natural rubber adhesives """
        return -70, "degC"
    @property
    def Tm(self):
        """ adhesives often degrade, no well-defined melt """
        return None
    def crystallinity(self, T=None):
        """Crystallinity of the solid phase"""
        return 0
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0


# -- AdhesiveSyntheticRubber ------------------------------------------
class AdhesiveSyntheticRubber(layer):
    _chemicalsubstance = "cis-1,4-polyisoprene" # styrene-butadiene rubber (SBR) or similar
    _polarityindex = 2.0  # non-polar or slightly polar, depending on rubber type.
    """ extended pantankar.layer for synthetic rubber adhesives """
    def __init__(self, l=20e-6, D=1e-13, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="adhesive synthetic rubber",**extra):
        """ constructor for a synthetic rubber-based adhesive layer """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="adhesive",
            layermaterial="synthetic rubber adhesive",
            layercode="sRubber",
            **extra
        )
    def density(self, T=None):
        """ typical density ~920 kg/m^3, adjust as needed """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 920 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"
    @property
    def Tg(self):
        """ approximate Tg of synthetic rubber adhesives """
        return -50, "degC"
    @property
    def Tm(self):
        """ adhesives often degrade, no well-defined melt """
        return None
    def crystallinity(self, T=None):
        """Crystallinity of the solid phase"""
        return 0
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0


# -- AdhesiveEVA (ethylene-vinyl acetate) ------------------------------
class AdhesiveEVA(layer):
    _chemicalsubstance = "ethylene" # Ethylene + vinyl acetate
    _polarityindex = 2.5  # Mostly non-polar backbone with some polar acetate groups.
    """ extended pantankar.layer for EVA-based adhesives """
    def __init__(self, l=20e-6, D=1e-13, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="adhesive EVA",**extra):
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="adhesive",
            layermaterial="EVA adhesive",
            layercode="EVA",
            **extra
        )
    def density(self, T=None):
        """ typical density ~930 kg/m^3 """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 930 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"
    @property
    def Tg(self):
        """ approximate Tg of EVA adhesives """
        return -30, "degC"
    @property
    def Tm(self):
        """ adhesives often degrade, no well-defined melt """
        return None

    def crystallinity(self, T=None):
        """Crystallinity of the solid phase"""
        return 0
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0


# -- AdhesiveVAE (vinyl acetate-ethylene) -----------------------------
class AdhesiveVAE(layer):
    _chemicalsubstance = "vinyl acetate" # Ethylene + vinyl acetate
    _polarityindex = 4.0  # More polar than EVA (larger fraction of acetate).
    """ extended pantankar.layer for VAE adhesives """
    def __init__(self, l=20e-6, D=1e-13, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="adhesive VAE",**extra):
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="adhesive",
            layermaterial="VAE adhesive",
            layercode="VAE",
            **extra
        )
    def density(self, T=None):
        """ typical density ~950 kg/m^3 """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 950 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"
    @property
    def Tg(self):
        """ approximate Tg of VAE adhesives """
        return 10, "degC"
    @property
    def Tm(self):
        """ adhesives often degrade, no well-defined melt """
        return None
    def crystallinity(self, T=None):
        """Crystallinity of the solid phase"""
        return 0
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0


# -- AdhesivePVAC (polyvinyl acetate) ---------------------------------

class AdhesivePVAC(layer):
    """ extended pantankar.layer for PVAc adhesives """
    _chemicalsubstance = "vinyl acetate" # Vinyl acetate (CH₂=CHO–Ac)
    _polarityindex = 7.0  # PVAc is fairly polar (acetate groups)
    def __init__(self, l=20e-6, D=1e-13, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="adhesive PVAc",**extra):
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="adhesive",
            layermaterial="PVAc adhesive",
            layercode="PVAc",
            **extra
        )
    def density(self, T=None):
        """ typical density ~1100 kg/m^3 """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 1100 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"
    @property
    def Tg(self):
        """ approximate Tg of PVAc adhesives """
        return 35, "degC"
    @property
    def Tm(self):
        """ adhesives often degrade, no well-defined melt """
        return None
    def crystallinity(self, T=None):
        """Crystallinity of the solid phase"""
        return 0
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0


# -- AdhesiveAcrylate -------------------------------------------------
class AdhesiveAcrylate(layer):
    """ extended pantankar.layer for acrylate adhesives """
    _chemicalsubstance = "n-butyl acrylate" # Acrylic esters (e.g. n-butyl acrylate)
    _polarityindex = 6.0  # Ester groups confer moderate polarity
    def __init__(self, l=20e-6, D=1e-13, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="adhesive acrylate",**extra):
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="adhesive",
            layermaterial="acrylate adhesive",
            layercode="Acryl",
            **extra
        )
    def density(self, T=None):
        """ typical density ~1000 kg/m^3 """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 1000 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"
    @property
    def Tg(self):
        """ approximate Tg of acrylate adhesives """
        return -20, "degC"
    @property
    def Tm(self):
        """ adhesives often degrade, no well-defined melt """
        return None
    def crystallinity(self, T=None):
        """Crystallinity of the solid phase"""
        return 0
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0


# -- AdhesivePU (polyurethane) ----------------------------------------
class AdhesivePU(layer):
    """ extended pantankar.layer for polyurethane adhesives """
    _chemicalsubstance = "diisocyanate" # Diisocyanate + polyol (–NH–CO–O–)
    _polarityindex = 5.0  # Can vary widely by chemistry; moderate polarity.
    def __init__(self, l=20e-6, D=1e-13, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="adhesive PU",**extra):
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="adhesive",
            layermaterial="polyurethane adhesive",
            layercode="PU",
            **extra
        )
    def density(self, T=None):
        """ typical density ~1100 kg/m^3 """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 1100 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"
    @property
    def Tg(self):
        """ approximate Tg of polyurethane adhesives """
        return -50, "degC"
    @property
    def Tm(self):
        """ adhesives often degrade, no well-defined melt """
        return None
    def crystallinity(self, T=None):
        """Crystallinity of the solid phase"""
        return 0
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0


# <<<<<<<<<<<<<<<<<<<<<<< P A P E R   &   C A R D B O A R D >>>>>>>>>>>>>>>>>>>>>>

# -- Paper ------------------------------------------------------------
class Paper(layer):
    """ extended pantankar.layer for paper (cellulose-based) """
    _physicalstate = "porous"       # solid (default), liquid, gas, porous
    _chemicalclass = "other"       # polymer (default), other
    _chemicalsubstance = "cellulose" # Cellulose (β-D-glucopyranose units)
    _polarityindex = 8.5  # Highly polar, strong hydrogen-bonding.
    def __init__(self, l=80e-6, D=1e-15, T=None,  # a guess for barrier properties
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="paper layer",**extra):
        """ Paper layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="paper",
            layermaterial="paper",
            layercode="paper",
            **extra
        )
    def density(self, T=None):
        """
        approximate density for typical paper ~800 kg/m^3
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 800 * (1 - 3*(T - layer._defaults["Td"]) * 1e-5), "kg/m**3"
    @property
    def Tg(self):
        """
        glass transition temperature is not typically used for paper,
        but we provide a placeholder.
        """
        return 200, "degC"  # purely illustrative placeholder
    @property
    def Tm(self):
        """ paper is not a pure polymer melt -> no Tm """
        return None
    def crystallinity(self, T=None):
        """Crystallinity of the solid phase"""
        return 1
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0.9

# -- Cardboard --------------------------------------------------------
class Cardboard(layer):
    """ extended pantankar.layer for cardboard (cellulose-based) """
    _physicalstate = "porous"       # solid (default), liquid, gas, porous
    _chemicalclass = "other"       # polymer (default), other
    _chemicalsubstance = "cellulose"
    _polarityindex = 8.0  # Can vary widely by chemistry; moderate polarity.
    def __init__(self, l=500e-6, D=1e-15, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="cardboard layer",**extra):
        """ Cardboard layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="paper",
            layermaterial="cardboard",
            layercode="board",
            **extra
        )
    def density(self, T=None):
        """
        approximate density for typical cardboard ~700 kg/m^3
        """
        T = self.T[0] if T is None else check_units(T, None, "degC")[0]
        return 700 * (1 - 3*(T - layer._defaults["Td"]) * 1e-5), "kg/m**3"
    @property
    def Tg(self):
        """
        same placeholder concept for paper-based material
        """
        return 200, "degC"
    @property
    def Tm(self):
        """ cardboard also no real melt """
        return None
    def crystallinity(self, T=None):
        """Crystallinity of the solid phase"""
        return 1
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 0.8

# <<<<<<<<<<<<<<<<<<<<<<< G A S E S  >>>>>>>>>>>>>>>>>>>>>>

# <-- air | ideal gas layer ---------------------------------->
class air(layer):
    """  extended pantankar.layer for ideal gases such as air """
    _physicalstate = "gas"       # solid (default), liquid, gas, porous
    _chemicalclass = "other"       # polymer (default), other
    def __init__(self,l=1e-2,D=1e-6,T=None,
                 lunit=None,Dunit=None,Cunit=None,
                 layername="air layer",layercode="air",**extra):
        """ air layer constructor """
        T = layer._defaults["T"] if T is None else check_units(T,None,"degC")[0]
        TK = constants["T0K"]+T
        kair = 1/(constants["R"] *TK)
        kairunit = constants["iRT0Kunit"]
        layer.__init__(self,
                       l=l,D=D,k=kair,C0=0,T=T,
                       lunit=lunit,Dunit=Dunit,kunit=kairunit,Cunit=Cunit,
                       layername=layername,
                       layertype="air", # set by default at inititialization
                       layermaterial="ideal gas",
                       layercode="gas",
                       **extra
                       )
    def density(self, T=None):
        """Density of air at atmospheric pressure: density(T in K)"""
        TK = self.TK if T is None else check_units(T,None,"K")[0]
        P_atm = 101325  # Pa (1 atm)
        M_air = 28.9647 # g/mol = 0.0289647 kg/mol (Molar mass of dry air).
        return P_atm / ((constants["R"]/M_air) * TK), "kg/m**3"
    @property
    def Tm(self):
        """ air does not have a melting point in this context """
        return None
    def crystallinity(self, T=None):
        """Crystallinity of the solid phase"""
        return 1
    @property
    def porosity(self):
        """Porosity: volume fraction of voids"""
        return 1


# %% For testing and debugging
# ===================================================
# main()
# ===================================================
# for debugging purposes (code called as a script)
# the code is called from here
# ===================================================
if __name__ == '__main__':
    from patankar.loadpubchem import migrant

    # debug
    A = layer()
    D = layerLink("D")
    D[0]=1.2345e-12
    A.D = D
    repr(A)



    list_of_migrants = ["toluene","limonene","BHT","Irganox 1076"]
    m = [migrant(s) for s in list_of_migrants]
    P = PS();
    P.update(migrant=m[1]).D
    P.update(migrant=m[2]).D


    materials = gPET()+wPET()+PS()+PP()
    m = [migrant(s) for s in list_of_migrants]
    materials.update(migrant=m[0],T=50)
    print(materials.D)
    materials.update(migrant=m[1])


    P=PS(migrant="toluene",T=20)
    P.D
    repr(P)


    P = rPS(migrant="bisphenol A",T=40)
    P.D
    repr(P)

    P=wPET(migrant="toluene",T=40)
    P.D
    repr(P)

    G = air(T=60)
    P = LDPE(D=1e-8,Dunit='cm**2/s')
    P = LDPE(D=(1e-8,"cm**2/s"))
    A = LDPE()
    A=layer(D=1e-14,l=50e-6)
    print("\n",repr(A),"\n"*2)
    A
    B=A*3
    D = B[1:2]
    B=A+A
    C=B[1]
    B.l = [1,2]
    A.struct()
    E = B*4
    #E[1:4]=[]
    E
    # ======
    A = layer(layername = "layer A")
    B = layer(layername = "layer B")
    C = layer(layername = "layer C")
    D = layer(layername = "layer D")
    # test = A+B+C+D
    # test[2] = test[0]
    # test[3] = []
    # test
    test = A+A+B+B+B+C
    print("\n",repr(test),"\n"*2)
    testsimple = test.simplify()
    print("\n",repr(testsimple),"\n"*2)
    testsimple.mesh()

    # test with substance
    m1 = migrant(name='limonene')
    m2 = migrant(name='anisole')
    pet_with_limonene = gPET(substance=m1,D=None,T=40,l=(50,"um"))
    PP_with_anisole = PP(substance=m2,D=None,T=40,l=(200,"um"))
    print("\n",repr(pet_with_limonene),"\n"*2)

    test = pet_with_limonene + PP_with_anisole
    test.D
    print("\n",repr(test),"\n"*2)
