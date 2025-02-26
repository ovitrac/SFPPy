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

List of implemnted materials:

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


@version: 1.2
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2022-02-21
@rev. 2025-02-14

"""

# ---- History ----
# Created on Tue Jan 18 09:14:34 2022
# 2022-01-19 RC
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


# %% Private functions

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




# %% Core class: layer
# default values (usable in layer object methods)
# these default values can be moved in a configuration file


# Main class definition
# =======================
class layer:
    """
        layer class from patankar package
        ...
        strings properties: layername, layertype, layermaterial
        scalar properties: D,k,l,C0
        ...
        Example:
            A = layer(D=1e-14,l=50e-6,layername="layer A",layertype="polymer",layermaterial="PP")
    """
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
                 substance = None, Dmodel = None, kmodel = None,
                 nmesh=None, nmeshmin=None,
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

        # set substance and property models
        self._substance = substance

        # set history for all layers merged with +
        self._layerclass_history = []

        # set verbosity attributes
        self.verbosity = 0 if verbosity is None else verbosity
        self.verbose = verbosity>0 if verbose is None else verbose

        # we initialize the acknowlegment process for future property propagation
        self._hasbeeninherited = {}

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
    # Class method to handle ambiguous definition from end-user
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
    def __add__(self,other):
        """ C=A+B | overload + operator """
        if isinstance(other, layer):
            res = duplicate(self)
            res._nmeshmin = min(self._nmeshmin,other._nmeshmin)
            # propage substance
            if self._substance is None:
                res._substance = other._substance
            else:
                if isinstance(self._substance,migrant) and isinstance(other._substance,migrant):
                    if self._substance.M != other._substance.M:
                        print("Warning: the smallest subtance is propagated everywhere")
                    res._substance = self._substance if self._substance.M<=other._substance.M else other._substance
                else:
                    res._substance = None
            for p in ["_name","_type","_material","_code","_nlayer"]:
                setattr(res,p,getattr(self,p)+getattr(other,p))
            for p in ["_l","_D","_k","_C0","_rho","_T"]:
                setattr(res,p,np.concatenate((getattr(self,p),getattr(other,p))))
            # we add the history of all layers
            res._layerclass_history = self.layerclass_history + other.layerclass_history
            return res
        else: raise ValueError("invalid layer object")


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
        for p in ["_name","_type","_material","_l","_D","_k","_C0"]:
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
            for p in ["_name","_type","_material","_l","_D","_k","_C0"]:
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
            for p in ["_name","_type","_material","_l","_D","_k","_C0"]:
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
    @property
    def layerclass_history(self):
        return self._layerclass_history if self._layerclass_history != [] else [self.layerclass]
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
    def l(self): return self._l
    @property
    def D(self):
        Dtmp = self.Dmodel()
        if Dtmp is not None:
            return np.full_like(self._D, Dtmp)
        else:
            return self._D
    @property
    def k(self):
        ktmp = self.kmodel()
        if ktmp:
            return np.full_like(self._D, ktmp)
        else:
            return self._k
    @property
    def C0(self): return self._C0
    @property
    def rho(self): return self._rho
    @property
    def T(self): return self._T
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
    def rank(self): return np.flip(np.argsort(np.array(self.resistance))+1).tolist()
    @property
    def referencelayer(self): return np.argmax(self.resistance)
    @property
    def lreferencelayer(self): return self.l[self.referencelayer]
    @property
    def Foscale(self): return self.D[self.referencelayer]/self.lreferencelayer**2

    # layer substance (of class migrant or None)
    @property
    def substance(self): return self._substance

    # Dmodel and kmodel returned as properties (they are lambda functions)
    # polymer and mass are udpdated on the fly (the code loops over all layers)
    @property
    def Dmodel(self):
        """Return a callable function that evaluates D with updated parameters."""
        if not isinstance(self._substance,migrant) or self._substance.Deval() is None:
            return lambda **kwargs: None  # Return a function that always returns None
        template = self._substance.Dtemplate.copy()
        template.update()
        def func(**kwargs):
            D = np.empty_like(self._D)
            for (i,),T in np.ndenumerate(self.T.ravel()): # loop over all layers via T
                template.update(polymer=self.layerclass_history[i],T=T) # updated layer properties
                # inherit eventual user parameters
                D[i] = self._substance.D.evaluate(**dict(template, **kwargs))
            return D
        return func

    @property
    def kmodel(self):
        """Return a callable function that evaluates k with updated parameters."""
        if not isinstance(self._substance,migrant) or self._substance.keval() is None:
            return lambda **kwargs: None  # Return a function that always returns None
        template = self._substance.ktemplate.copy()
        template.update()
        def func(**kwargs):
            k = np.empty_like(self._k)
            for (i,),T in np.ndenumerate(self.T.ravel()): # loop over all layers via T
                template.update(polymer=self.layerclass_history[i],T=T) # updated layer properties
                # inherit eventual user parameters
                k[i] = self._substance.k.evaluate(**dict(template, **kwargs))
            return k
        return func


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
    # Getter methods and tools to validate inputs checknumvalue and checktextvalue
    # --------------------------------------------------------------------
    def checknumvalue(self,value,ExpectedUnits=None):
        """ returns a validate value to set properties """
        if isinstance(value,tuple):
            value = check_units(value,ExpectedUnits=ExpectedUnits)
        if isinstance(value,int): value = float(value)
        if isinstance(value,float): value = np.array([value])
        if isinstance(value,list): value = np.array(value)
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
        if not isinstance(value,migrant):
            raise TypeError(f"value must be a migrant class not a {type(value).__name__}")
        self._substance = value


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
            hasDmodel = self.Dmodel() is not None
            haskmodel = self.kmodel() is not None
            properties_ = {"l":False,"D":hasDmodel,"k":haskmodel,"C0":False}
            if hasDmodel or haskmodel:
                properties_["T"] = False
            fmtval = '%10s: '+self._printformat+" [%s]"
            fmtstr = '%10s= %s'
            if self._nlayer==1:
                print(f'monolayer of {self.__description}:')
            else:
                print(f'{self._nlayer}-multilayer of {self.__description}:')
            for n in range(1,self._nlayer+1):
                modelinfo = {
                    "D": f"{self._substance.D.__name__}({self.layerclass_history[n-1]},{self._substance},T={float(self.T[0])} {self.Tunit})" if hasDmodel else "",
                    "k": f"{self._substance.k.__name__}({self.layerclass_history[n-1]},{self._substance} g/mol,T={float(self.T[0])} {self.Tunit})" if haskmodel else "",
                    }
                print('-- [ layer %d of %d ] ---------- barrier rank=%d --------------'
                      % (n,self._nlayer,self.rank[n-1]))
                for p in ["name","type","material","code"]:
                    v = getattr(self,p)
                    print('%10s: "%s"' % (p,v[n-1]),flush=True)
                for p in properties_.keys():
                    v = getattr(self,p)                 # value
                    vunit = getattr(self,p[0]+"unit")   # value unit
                    print(fmtval % (p,v[n-1],vunit),flush=True)
                    if properties_[p]:
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
    # --------------------------------------------------------------------
    def _from(self,medium=None):
        """Propagates contact conditions from food instance"""
        from patankar.food import foodphysics
        if not isinstance(medium,foodphysics):
            raise TypeError(f"medium must be a foodphysics, foodlayer not a {type(medium).__name__}")
        if not hasattr(medium, "contacttemperature"):
            medium.contacttemperature = self.T[0]
        T = self.get_param("contacttemperature",40,acceptNone=False)
        self.T = np.full_like(self.T,T)

    # overloading operation
    def __lshift__(self, medium):
        """Overloads << to propagate contact conditions from food."""
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
        return self.migration(medium,**kwargs)

    def migration(self,medium,**kwargs):
        from patankar.migration import senspatankar
        sim = senspatankar(self,medium,**kwargs)
        medium.lastsimulation = sim # store the last simulation result in medium
        medium.lastinput = self # store the last input (self)
        sim.savestate(self,medium) # store store the inputs in sim for chaining
        return sim

    # overloading operation
    def __rshift__(self, medium):
        """Overloads >> to propagate migration to food."""
        from patankar.food import foodphysics
        if not isinstance(medium,foodphysics):
            raise TypeError(f"medium must be a foodphysics object not a {type(medium).__name__}")
        return self.contact(medium)

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

=========================================================
"""

# <<<<<<<<<<<<<<<<<<<<<<< P O L Y O L E F I N S >>>>>>>>>>>>>>>>>>>>>>

# <-- LDPE polymer ---------------------------------->
class LDPE(layer):
    """  extended pantankar.layer for low-density polyethylene LDPE  """
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
        T = self.T if T is None else check_units(T,None,"degC")[0]
        return 920 *(1-3*(T-layer._defaults["Td"])*20e-5),"kg/m**3" # lowest temperature
    @property
    def Tg(self):
        """ glass transition temperature of LDPE """
        return -130,"degC" # lowest temperature


# <-- HDPE polymer ---------------------------------->
class HDPE(layer):
    """  extended pantankar.layer for high-density polyethylene HDPE  """
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
        T = self.T if T is None else check_units(T,None,"degC")[0]
        return 940 *(1-3*(T-layer._defaults["Td"])*11e-5),"kg/m**3" # lowest temperature
    @property
    def Tg(self):
        """ glass transition temperature of HDPE """
        return -100,"degC" # highest temperature

# <-- LLDPE polymer ---------------------------------->
class LLDPE(layer):
    """ extended pantankar.layer for linear low-density polyethylene LLDPE """
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        # Similar formula to LDPE and HDPE, with a coefficient suitable for LLDPE.
        return 915 * (1 - 3 * (T - layer._defaults["Td"]) * 15e-5), "kg/m**3"
    @property
    def Tg(self):
        """
        glass transition temperature of LLDPE
        Typically close to LDPE, though slightly higher or lower can be found in the literature.
        """
        return -120, "degC"

# <-- PP polymer ---------------------------------->
class PP(layer):
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
        T = self.T if T is None else check_units(T,None,"degC")[0]
        return 910 *(1-3*(T-layer._defaults["Td"])*7e-5),"kg/m**3" # lowest temperature
    @property
    def Tg(self):
        """ glass transition temperature of PP """
        return 0,"degC" # highest temperature

# -- PPrubber (atactic polypropylene) ---------------------------------
class PPrubber(layer):
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 900 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"

    @property
    def Tg(self):
        """ glass transition temperature of atactic/rubbery PP """
        return -20, "degC"


# -- oPP (bioriented polypropylene) ------------------------------------
class oPP(layer):
    """ extended pantankar.layer for bioriented polypropylene oPP """
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 910 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"

    @property
    def Tg(self):
        """ glass transition temperature of bioriented PP """
        return 0, "degC"


# <<<<<<<<<<<<<<<<<<<<<<< P O L Y V I N Y L S >>>>>>>>>>>>>>>>>>>>>>

# -- PS (polystyrene) -----------------------------------------------
class PS(layer):
    """ extended pantankar.layer for polystyrene (PS) """
    def __init__(self, l=100e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in PS",**extra):
        """ PS layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="polymer",
            layermaterial="polystyrene",
            layercode="PS",
            **extra
        )
    def density(self, T=None):
        """
        density of PS: ~1050 kg/m^3
        """
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 1050 * (1 - 3*(T - layer._defaults["Td"]) * 5e-5), "kg/m**3"

    @property
    def Tg(self):
        """ glass transition temperature of PS """
        return 100, "degC"


# -- HIPS (high-impact polystyrene) -----------------------------------
class HIPS(layer):
    """ extended pantankar.layer for high-impact polystyrene (HIPS) """
    def __init__(self, l=100e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in HIPS",**extra):
        """ HIPS layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="polymer",
            layermaterial="high-impact polystyrene",
            layercode="HIPS",
            **extra
        )
    def density(self, T=None):
        """
        density of HIPS: ~1040 kg/m^3
        """
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 1040 * (1 - 3*(T - layer._defaults["Td"]) * 5e-5), "kg/m**3"

    @property
    def Tg(self):
        """ glass transition temperature of HIPS """
        return 95, "degC"


# -- PBS (assuming a styrene-based polymer) ---------------------------
class SBS(layer):
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 1030 * (1 - 3*(T - layer._defaults["Td"]) * 5e-5), "kg/m**3"

    @property
    def Tg(self):
        """ glass transition temperature of 'DBS' """
        return 90, "degC"


# -- rigidPVC ---------------------------------------------------------
class rigidPVC(layer):
    """ extended pantankar.layer for rigid PVC """
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 1400 * (1 - 3*(T - layer._defaults["Td"]) * 5e-5), "kg/m**3"

    @property
    def Tg(self):
        """ glass transition temperature of rigid PVC """
        return 80, "degC"


# -- plasticizedPVC ---------------------------------------------------
class plasticizedPVC(layer):
    """ extended pantankar.layer for plasticized PVC """
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 1300 * (1 - 3*(T - layer._defaults["Td"]) * 5e-5), "kg/m**3"

    @property
    def Tg(self):
        """ glass transition temperature of plasticized PVC """
        return -40, "degC"


# <<<<<<<<<<<<<<<<<<<<<<< P O L Y E S T E R S >>>>>>>>>>>>>>>>>>>>>>

# -- gPET (glassy PET, T < 76°C) --------------------------------------
class gPET(layer):
    """ extended pantankar.layer for PET in its glassy state (below ~76°C) """
    def __init__(self, l=200e-6, D=1e-14, T=None,
                 k=None, C0=None, lunit=None, Dunit=None, kunit=None, Cunit=None,
                 layername="layer in gPET",**extra):
        """ glassy PET layer constructor """
        super().__init__(
            l=l, D=D, k=k, C0=C0, T=T,
            lunit=lunit, Dunit=Dunit, kunit=kunit, Cunit=Cunit,
            layername=layername,
            layertype="polymer",
            layermaterial="glassy PET",
            layercode="PET",
            **extra
        )
    def density(self, T=None):
        """
        density of glassy PET: ~1350 kg/m^3
        """
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 1350 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"

    @property
    def Tg(self):
        """ approximate glass transition temperature of PET """
        return 76, "degC"


# -- rPET (rubbery PET, T > 76°C) --------------------------------------
class rPET(layer):
    """ extended pantankar.layer for PET in its rubbery state (above ~76°C) """
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 1350 * (1 - 3*(T - layer._defaults["Td"]) * 1e-4), "kg/m**3"

    @property
    def Tg(self):
        """ approximate glass transition temperature of PET """
        return 76, "degC"


# -- PBT --------------------------------------------------------------
class PBT(layer):
    """ extended pantankar.layer for polybutylene terephthalate (PBT) """
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 1310 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"

    @property
    def Tg(self):
        """ glass transition temperature of PBT """
        return 40, "degC"


# -- PEN --------------------------------------------------------------
class PEN(layer):
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 1330 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"

    @property
    def Tg(self):
        """ glass transition temperature of PEN """
        return 120, "degC"


# <<<<<<<<<<<<<<<<<<<<<<< P O L Y A M I D E S >>>>>>>>>>>>>>>>>>>>>>

# -- PA6 --------------------------------------------------------------
class PA6(layer):
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 1140 * (1 - 3*(T - layer._defaults["Td"]) * 5e-5), "kg/m**3"

    @property
    def Tg(self):
        """ glass transition temperature of PA6 """
        return 50, "degC"


# -- PA66 -------------------------------------------------------------
class PA66(layer):
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 1150 * (1 - 3*(T - layer._defaults["Td"]) * 5e-5), "kg/m**3"

    @property
    def Tg(self):
        """ glass transition temperature of PA66 """
        return 70, "degC"


# <<<<<<<<<<<<<<<<<<<<<<< A D H E S I V E S >>>>>>>>>>>>>>>>>>>>>>

# -- AdhesiveNaturalRubber --------------------------------------------
class AdhesiveNaturalRubber(layer):
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 910 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"

    @property
    def Tg(self):
        """ approximate Tg of natural rubber adhesives """
        return -70, "degC"


# -- AdhesiveSyntheticRubber ------------------------------------------
class AdhesiveSyntheticRubber(layer):
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 920 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"

    @property
    def Tg(self):
        """ approximate Tg of synthetic rubber adhesives """
        return -50, "degC"


# -- AdhesiveEVA (ethylene-vinyl acetate) ------------------------------
class AdhesiveEVA(layer):
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 930 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"
    @property
    def Tg(self):
        """ approximate Tg of EVA adhesives """
        return -30, "degC"


# -- AdhesiveVAE (vinyl acetate-ethylene) -----------------------------
class AdhesiveVAE(layer):
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 950 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"
    @property
    def Tg(self):
        """ approximate Tg of VAE adhesives """
        return 10, "degC"


# -- AdhesivePVAC (polyvinyl acetate) ---------------------------------
class AdhesivePVAC(layer):
    """ extended pantankar.layer for PVAc adhesives """
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 1100 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"
    @property
    def Tg(self):
        """ approximate Tg of PVAc adhesives """
        return 35, "degC"


# -- AdhesiveAcrylate -------------------------------------------------
class AdhesiveAcrylate(layer):
    """ extended pantankar.layer for acrylate adhesives """
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 1000 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"
    @property
    def Tg(self):
        """ approximate Tg of acrylate adhesives """
        return -20, "degC"


# -- AdhesivePU (polyurethane) ----------------------------------------
class AdhesivePU(layer):
    """ extended pantankar.layer for polyurethane adhesives """
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 1100 * (1 - 3*(T - layer._defaults["Td"]) * 7e-5), "kg/m**3"
    @property
    def Tg(self):
        """ approximate Tg of polyurethane adhesives """
        return -50, "degC"


# <<<<<<<<<<<<<<<<<<<<<<< P A P E R   &   C A R D B O A R D >>>>>>>>>>>>>>>>>>>>>>

# -- Paper ------------------------------------------------------------
class Paper(layer):
    """ extended pantankar.layer for paper (cellulose-based) """
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 800 * (1 - 3*(T - layer._defaults["Td"]) * 1e-5), "kg/m**3"

    @property
    def Tg(self):
        """
        glass transition temperature is not typically used for paper,
        but we provide a placeholder.
        """
        return 200, "degC"  # purely illustrative placeholder


# -- Cardboard --------------------------------------------------------
class Cardboard(layer):
    """ extended pantankar.layer for cardboard (cellulose-based) """
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
        T = self.T if T is None else check_units(T, None, "degC")[0]
        return 700 * (1 - 3*(T - layer._defaults["Td"]) * 1e-5), "kg/m**3"

    @property
    def Tg(self):
        """
        same placeholder concept for paper-based material
        """
        return 200, "degC"





# <<<<<<<<<<<<<<<<<<<<<<< G A S E S  >>>>>>>>>>>>>>>>>>>>>>

# <-- air | ideal gas layer ---------------------------------->
class air(layer):
    """  extended pantankar.layer for ideal gases such as air """
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

# %% For testing and debugging
# ===================================================
# main()
# ===================================================
# for debugging purposes (code called as a script)
# the code is called from here
# ===================================================
if __name__ == '__main__':
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
