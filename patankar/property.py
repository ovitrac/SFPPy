#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
SFPPy Module: Property
===============================================================================
Defines physical parameters for mass transfer, independent of specific theoretical or empirical models.
Currently implements the **Piringer model** for worst-case migration simulations.

**Main Components:**
- **Base Class: `migrationProperty`** (Holds generic attributes for any mass transfer property)
- **Subclasses for Specific Properties:**
    - `Diffusivities`: Defines diffusion coefficients (D)
    - `HenriLikeCoeffcicients`: Defines Henry-like coefficients (k)
    - `ActivityCoeffcicients`: Defines activity coefficients (γ)
    - `PartitionCoeffcicients`: Defines partition coefficients (K)
- **Piringer Model (`Dpiringer`)**
    - Empirical overestimation model for polymer diffusion
    - Used for migration simulations in `migration.py`
    - Directly invoked by `loadpubchem.py` when retrieving substance properties

**Integration with SFPPy Modules:**
- Used by `loadpubchem.py` to predict missing diffusivity or partitioning values for chemical migrants.
- Applied in `migration.py` for solving mass transfer equations.

Example:
```python
from property patankar.import Dpiringer
D_value = Dpiringer.evaluate(polymer="LDPE", M=100, T=40)
```


===============================================================================
Details
===============================================================================
This module offers the necessary abstraction to any physical parameter governing
mass transfer independently of the applied molecular theory or emprical model used
to calculate them.

Currently this module implements seamlesly the Dpringer model for worst-case simulations
for risk assessment.

This class is used directly by loadpubchem without involving any user operation.
The name or CAS of the substance will trigger the predictions for the considered application.


@version: 1.2
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2024-07-21
@rev: 2025-02-21

"""

import math

# %% Top classes for any property
# level 0
class migrationProperty:
    """Base class to hold general properties used for migration of substances."""
    property = "any"
    notation = ""
    description = "root class"
    name = "root"
    parameters = []  # e.g. ["M", "T"]
    SIunits = ""

    # private properties
    _model = ""
    _theory = ""
    _source = ""
    _author = "olivier.vitrac@agroparistech.fr"
    _license = "MIT"
    _version = 0.1
    _available_to_import = False


    def __repr__(self):
        """Formatted string representation for nice display."""
        # Define attribute names and their corresponding values
        attributes = {
            "property": self.property,
            "notation": self.notation,
            "description": self.description,
            "name": self.name,
            "parameters": self.parameters,
            "SIunits": self.SIunits,
            "model": self._model,
            "theory": self._theory,
            "source": self._source,
            "author": self._author,
            "license": self._license,
            "version": self._version,
        }
        # Filter out None or empty string values
        filtered_attributes = {k: v for k, v in attributes.items() if v not in (None, "")}
        # Find the max length of attribute names for alignment
        max_key_length = max(len(k) for k in filtered_attributes.keys()) if filtered_attributes else 0
        # Format the output with proper alignment
        lines = [f"{k.rjust(max_key_length)}: {v}" for k, v in filtered_attributes.items()]
        print("\n".join(lines))
        return str(self)

    def __str__(self):
        """Formatted string representation of property"""
        return f"<{self.__class__.__name__}: {self.property}:{self.notation}>"

# level 1
class Diffusivities(migrationProperty):
    """Base class for diffusion-related models."""
    property = "Diffusivity"
    notation = "D"
    description = "Mathematical model to estimate diffusivities"
    SIunits = "m**2/s"

class HenriLikeCoeffcicients(migrationProperty):
    property = "Henri-like coefficient"
    notation = "k"
    description = "Mathematical model to estimate Henri-like coefficients"
    SIunits = None

class ActivityCoeffcicients(migrationProperty):
    property = "Activity coefficient"
    notation = "g"
    description = "Mathematical model to estimate activity coefficients"
    SIunits = None

class PartitionCoeffcicients(migrationProperty):
    property = "Partition Coefficient"
    notation = "K"
    description = "Mathematical model to estimate partition coefficients"
    SIunits = None

# %% Piringer model
class Dpiringer(Diffusivities):
    """
        Piringer's overestimate of diffusion coefficient.

        Two implementations are offered in the class:
            - static: Dpiringer.evaluate(polymer="polymer",M=Mvalue,T=Tvalue)
            - dynamic: Dmodel = Dpiringer(polymer="polymer"...)
                       Dmodel.eval(M=Mvalue,T=Tvalue)

    """
    name = "Piringer"
    description = "Piringer's overestimate of diffusion coefficients"
    model = "empirical"
    theory = "scaling"
    parameters = {"M": {"description": "molecular mass","units": "g/mol"},
                  "T": {"description": "temperature","units": "degC"}
                }
    _available_to_import = True # this model can be directly imported

    # Piringer values (the primary key matches the one used in layer)
    piringer_data = {
        # -- polyolefins -------------------------------------------
        "HDPE": {    # category: polyolefins
            "className": "HDPE",
            "type": "polymer",
            "material": "high-density polyethylene",
            "code": "HDPE",
            "description": "Piringer parameters for HDPE.",
            "App": 14.5,
            "tau": 1577
        },
        "LDPE": {    # category: polyolefins
            "className": "LDPE",
            "type": "polymer",
            "material": "low-density polyethylene",
            "code": "LDPE",
            "description": "Piringer parameters for LDPE.",
            "App": 11.5,
            "tau": 0
        },
        "LLDPE": {   # category: polyolefins
            "className": "LLDPE",
            "type": "polymer",
            "material": "linear low-density polyethylene",
            "code": "LLDPE",
            "description": "Piringer parameters for LLDPE.",
            "App": 11.5,
            "tau": 0
        },
        "PP": {      # category: polyolefins
            "className": "PP",
            "type": "polymer",
            "material": "isotactic polypropylene",
            "code": "PP",
            "description": "Piringer parameters for isotactic PP.",
            "App": 13.1,
            "tau": 1577
        },
        "aPP": {     # category: polyolefins
            "className": "PPrubber",
            "type": "polymer",
            "material": "atactic polypropylene",
            "code": "aPP",
            "description": "Piringer parameters for atactic PP.",
            "App": 11.5,
            "tau": 0
        },
        "oPP": {     # category: polyolefins
            "className": "oPP",
            "type": "polymer",
            "material": "bioriented polypropylene",
            "code": "oPP",
            "description": "Piringer parameters for bioriented PP.",
            "App": 13.1,
            "tau": 1577
        },

        # -- polyvinyls --------------------------------------------
        "pPVC": {    # category: polyvinyls
            "className": "plasticizedPVC",
            "type": "polymer",
            "material": "plasticized PVC",
            "code": "pPVC",
            "description": "Piringer parameters for plasticized PVC.",
            "App": 14.6,
            "tau": 0
        },
        "PVC": {     # category: polyvinyls
            "className": "rigidPVC",
            "type": "polymer",
            "material": "rigid PVC",
            "code": "PVC",
            "description": "Piringer parameters for rigid PVC.",
            "App": -1.0,
            "tau": 0
        },

        # -- polystyrene, etc. (misc) ------------------------------
        "HIPS": {    # category: polystyrenics
            "className": "HIPS",
            "type": "polymer",
            "material": "high-impact polystyrene",
            "code": "HIPS",
            "description": "Piringer parameters for HIPS.",
            "App": 1.0,
            "tau": 0
        },
        "PBS": {     # category: polystyrenics
            "className": "PBS",
            "type": "polymer",
            "material": "styrene-based polymer PBS",
            "code": "PBS",
            "description": "No original Piringer data; set to None.",
            "App": 10.5,
            "tau": 0
        },
        "PS": {      # category: polystyrenics
            "className": "PS",
            "type": "polymer",
            "material": "polystyrene",
            "code": "PS",
            "description": "Piringer parameters for PS.",
            "App": -1.0,
            "tau": 0
        },

        # -- polyesters --------------------------------------------
        "PBT": {     # category: polyesters
            "className": "PBT",
            "type": "polymer",
            "material": "polybutylene terephthalate",
            "code": "PBT",
            "description": "Piringer parameters for PBT.",
            "App": 6.5,
            "tau": 1577
        },
        "PEN": {     # category: polyesters
            "className": "PEN",
            "type": "polymer",
            "material": "polyethylene naphthalate",
            "code": "PEN",
            "description": "Piringer parameters for PEN.",
            "App": 5.0,
            "tau": 1577
        },
        "PET": {     # category: polyesters
            "className": "gPET",
            "type": "polymer",
            "material": "glassy PET",
            "code": "PET",
            "description": "Piringer parameters for glassy PET (inf Tg).",
            "App": 3.1,
            "tau": 1577
        },
        "rPET": {    # category: polyesters
            "className": "rPET",
            "type": "polymer",
            "material": "rubbery PET",
            "code": "rPET",
            "description": "Piringer parameters for rubbery PET (sup Tg).",
            "App": 6.4,
            "tau": 1577
        },

        # -- polyamides --------------------------------------------
        "PA6": {     # category: polyamides
            "className": "PA6",
            "type": "polymer",
            "material": "polyamide 6",
            "code": "PA6",
            "description": "Piringer parameters for polyamide 6.",
            "App": 0.0,
            "tau": 0
        },
        "PA6,6": {   # category: polyamides
            "className": "PA66",
            "type": "polymer",
            "material": "polyamide 6,6",
            "code": "PA6,6",
            "description": "Piringer parameters for polyamide 6,6.",
            "App": 2.0,
            "tau": 0
        },

        # -- adhesives --------------------------------------------
        "Acryl": {  # category: adhesives
            "className": "AdhesiveAcrylate",
            "type": "adhesive",
            "material": "acrylate adhesive",
            "code": "Acryl",
            "description": "Piringer parameters for acrylate adhesive.",
            "App": 4.5,
            "tau": 83
        },
        "EVA": {    # category: adhesives
            "className": "AdhesiveEVA",
            "type": "adhesive",
            "material": "EVA adhesive",
            "code": "EVA",
            "description": "Piringer parameters for EVA adhesive.",
            "App": 6.6,
            "tau": -1270
        },
        "rubber": { # category: adhesives
            "className": "AdhesiveNaturalRubber",
            "type": "adhesive",
            "material": "natural rubber adhesive",
            "code": "rubber",
            "description": "Piringer parameters for natural rubber adhesive.",
            "App": 11.3,
            "tau": -421
        },
        "PU": {     # category: adhesives
            "className": "AdhesivePU",
            "type": "adhesive",
            "material": "polyurethane adhesive",
            "code": "PU",
            "description": "Piringer parameters for polyurethane adhesive.",
            "App": 4.0,
            "tau": 250
        },
        "PVAc": {   # category: adhesives
            "className": "AdhesivePVAC",
            "type": "adhesive",
            "material": "PVAc adhesive",
            "code": "PVAc",
            "description": "Piringer parameters for PVAc adhesive.",
            "App": 6.6,
            "tau": -1270
        },
        "sRubber": {  # category: adhesives
            "className": "AdhesiveSyntheticRubber",
            "type": "adhesive",
            "material": "synthetic rubber adhesive",
            "code": "sRubber",
            "description": "Piringer parameters for synthetic rubber adhesive.",
            "App": 11.3,
            "tau": -421
        },
        "VAE": {      # category: adhesives
            "className": "AdhesiveVAE",
            "type": "adhesive",
            "material": "VAE adhesive",
            "code": "VAE",
            "description": "Piringer parameters for VAE adhesive.",
            "App": 6.6,
            "tau": -1270
        },

        # -- paper and board ---------------------------------------
        "board_polar": {   # category: paper_and_board
            "className": "Cardboard",
            "type": "paper",
            "material": "cardboard",
            "code": "board",
            "description": "Piringer parameters for cardboard (polar migrants).",
            "App": 4,
            "tau": -1511
        },

        "board_apol": {   # category: paper_and_board
            "className": "Cardboard",
            "type": "paper",
            "material": "cardboard",
            "code": "board",
            "description": "Piringer parameters for cardboard (variant for apolar).",
            "App": 7.4,
            "tau": -1511
        },

        "paper": {   # category: paper_and_board
            "className": "Paper",
            "type": "paper",
            "material": "paper",
            "code": "paper",
            "description": "Piringer parameters for paper.",
            "App": 6.6,
            "tau": -1900
        },

        # -- air ----------------------------------------------------
        "gas": {     # category: air
            "className": "air",
            "type": "air",
            "material": "ideal gas",
            "code": "gas",
            "description": "No Piringer data for air; set to None.",
            "App": None,
            "tau": None
        }
    }

    def __init__(self, polymer="LDPE", M=100, T=40):
        """
        Instantiate a Dpiringer object for a specific polymer key
        (e.g. 'LDPE', 'PET', or 'air'). The corresponding App and tau
        are looked up and stored as instance attributes.
        """
        polymer_str = polymer.strip()
        if polymer_str not in self.piringer_data:
            print(f"No exact match for polymer key: {polymer_str!r}")

        params = Dpiringer.get_piringer_params(polymer_str)
        if params["App"] is None or params["tau"] is None:
            raise ValueError(f"Piringer parameters not defined (App or tau is None) for {polymer_str!r}")
        self._polymer = polymer_str
        self._M = M
        self._T = T
        self._App = params["App"]
        self._tau = params["tau"]

    @property
    def polymer(self) -> str:
        """Return the stored polymer code (e.g. 'PET')."""
        return self._polymer

    @property
    def App(self) -> float:
        """Piringer's App constant for the selected polymer."""
        return self._App

    @property
    def tau(self) -> float:
        """Piringer's tau constant for the selected polymer."""
        return self._tau

    @property
    def M(self) -> float:
        """Molecular mass of the solute."""
        return self._M

    @property
    def T(self) -> float:
        """Temperature in degC."""
        return self._T

    @M.setter
    def M(self,value): self._M = value

    @T.setter
    def T(self,value): self._T = value

    def eval(self, M=None, T=None, **extra):
        """
        Compute Piringer D for this polymer (already stored in the instance)
        at molecular mass M (g/mol) and temperature T (°C).
        """
        M = self._M if M is None else M
        T = self._T if T is None else T
        # Convert T (°C) to T (K)
        TK = T + 273.15
        # Piringer expression for D in m^2/s
        exponent = (self._App
                    - (self._tau / TK)
                    - 0.135 * (M ** (2.0 / 3.0))
                    + 0.003 * M
                    - 10454.0 / TK)
        return math.exp(exponent)


    @classmethod
    def get_piringer_params(cls,polymer: str, data: dict = piringer_data):
        """
        Look up an entry in piringer_data by:
          1) Dictionary key (e.g. "LDPE")
          2) 'code' field (e.g. "LDPE")
          3) 'className' field (e.g. "LDPE")

        The matching is done case-insensitively.

        - If an exact match is found in either of those fields, return that entry.
        - If no exact match is found, attempt partial matches across all three fields and
          display them in a neat Markdown table if multiple partial matches appear.
        - If none found or the data is incomplete (App or tau is None), raise ValueError.
        """
        if polymer is None:
            raise ValueError("Please provide a polymer/material name")
        if not isinstance(polymer,str):
            raise TypeError(f"polymer must be a str not a {type(polymer).__name__}")
        polymer_str = polymer.strip().lower()

        # STEP 1: Try to find a single exact match
        #         across (dict key) or (entry["code"]) or (entry["className"])
        matched_key = None
        for k, info in data.items():
            # Check dictionary key, code, className
            if (
                polymer_str == k.lower() or
                (info["code"] and polymer_str == info["code"].lower()) or
                (info["className"] and polymer_str == info["className"].lower())
            ):
                matched_key = k
                break

        if matched_key is not None:
            # We found an exact match.
            entry = data[matched_key]
            return entry

        # STEP 2: No exact match => build partial match candidates
        partial_matches = []
        for k, info in data.items():
            k_l = k.lower()
            c_l = info["code"].lower() if info["code"] else ""
            n_l = info["className"].lower() if info["className"] else ""
            if (polymer_str in k_l) or (polymer_str in c_l) or (polymer_str in n_l):
                partial_matches.append(k)

        if not partial_matches:
            # No partial matches
            raise ValueError(f"No match or suggestion found for '{polymer}'.")

        if len(partial_matches) == 1:
            # Only one partial match => treat it like an exact match
            matched_key = partial_matches[0]
            entry = data[matched_key]
            return entry

        # STEP 3: Multiple partial matches => show a table
        # We'll build a dynamic Markdown table with columns:
        #   Key | className | code | material
        suggestions = []
        for pm in partial_matches:
            info = data[pm]
            suggestions.append([
                pm, info["className"], info["code"], info["material"]
            ])

        # Headers
        headers = ["Key", "className", "code", "material"]
        # Find maximum width for each column
        col_widths = [len(h) for h in headers]
        for row in suggestions:
            for i, cell in enumerate(row):
                col_widths[i] = max(col_widths[i], len(cell))

        # Build the header row
        header_line = "| " + " | ".join(
            headers[i].ljust(col_widths[i]) for i in range(len(headers))
        ) + " |"

        # Separator
        sep_line = "|-" + "-|-".join("-" * w for w in col_widths) + "-|"

        # Rows
        row_lines = []
        for row in suggestions:
            row_line = "| " + " | ".join(
                row[i].ljust(col_widths[i]) for i in range(len(row))
            ) + " |"
            row_lines.append(row_line)

        markdown_table = "\n".join([header_line, sep_line] + row_lines)

        raise ValueError(
            f"No exact match found for '{polymer}'. "
            f"Possible partial matches:\n\n{markdown_table}"
        )

    # static method (alternative for one shot evaluation)
    @classmethod
    def evaluate(cls, polymer="LLDPE", M=100.0, T=40.0, **extra):
        """
        Evaluate D (Piringer) for a single polymer, molecular mass (M), and temperature (T in °C).
        Replicates the essential logic of the original MATLAB Dpiringer function.
        No vectorization is performed (handles one polymer at a time).

        Parameters
        ----------
        polymer : str
            Polymer name (e.g. 'LLDPE', 'LDPE', 'rPET', etc.) as listed in the original data structure.
        M : float
            Molecular mass (g/mol), default = 100.
        T : float
            Temperature in °C, default = 40.

        Returns
        -------
        float
            The estimated diffusion coefficient in m^2/s (Piringer's overestimate).
        """

        # get Piringer model parameter
        params = cls.get_piringer_params(polymer)
        if params["App"] is None or params["tau"] is None:
            raise ValueError(
                f"Data for '{polymer}' is incomplete: App or tau is None."
            )

        # Convert T (°C) to T (K)
        TK = T + 273.15

        # Compute Ap = App - tau / TK
        App = params['App']
        tau = params['tau']
        Ap  = App - tau / TK  # dimensionless exponent part

        # Piringer expression for D in m^2/s
        # D = exp( Ap - 0.135 * M^(2/3) + 0.003 * M - 10454 / TK )
        exponent = Ap - 0.135 * (M ** (2.0 / 3.0)) + 0.003 * M - 10454.0 / TK
        D = math.exp(exponent)
        return D

# %% Available models to expose to layer.py and food.py
# List below the importable models (currently only D, k, K models are possible)
# They are imported via:
#    from property import MigrationPropertyModels, MigrationPropertyModel_validator
# Note that the name of the attribute (eg, "Piringer") must match class.name (eg, Dpiringer.name)
# A strict validator is proposed as MigrationPropertyModel_validator()

MigrationPropertyModels = {
    "D":{
        "Piringer": Dpiringer
        },
    "k":{
        },
    "g":{
        },
    "K":{
        },
    }

# Function helper to get a strict control on property models used by layer.py and food.py
def MigrationPropertyModel_validator(model=None,name=None,notation=None):
    """ Returns True if the proposed model is valid for the requested migraton property """
    rootclass = "migrationProperty"
    expectedpropclass = {"D":"Diffusivities",
                     "k":"HenriLikeCoefficients",
                     "g":"ActivityCoefficients",
                     "K":"PartitionCoefficients"}

    def get_root_parent(cls,level):
        """Returns the root parent class just after 'object'."""
        mro = cls.mro()  # Get the Method Resolution Order (MRO)
        for base in mro[level:]:  # Skip the class itself
            if base is not object:
                return base.__name__
        return None  # If no valid parent found

    if model is None or name is None or notation is None:
        raise ValueError("model, name and notation are mandatory.")
    if notation not in MigrationPropertyModels:
        raise ValueError(f"the property {notation} is not defined in MigrationPropertyModels")
    if type(model).__name__!="type":
        raise TypeError(f"model should be a class (e.g., Dpiringer) not a {type(model).__name__}")
    if get_root_parent(model,2)!=rootclass:
        raise TypeError(f'model "{model.__name__}" is not of class migrationProperty')
    if get_root_parent(model,1)!=expectedpropclass[notation]:
        raise TypeError(f'model "{model.__name__}" is not of class {expectedpropclass[notation]}, but of class {get_root_parent(model,1)}')
    if not model._available_to_import:
        raise TypeError(f'model "{model.__name__}" is not flagged for import')
    if model.name!=name:
        raise ValueError(f'model name "{model.name}" does not match the supplied name "{name}"')
    if model.notation!=notation:
        raise ValueError(f'model notation "{model.notation}" does not match the supplied name "{notation}"')
    return True # if all tests passed


# %% debug and standalone
# -------------------------------------------------------------------
# Example usage (for debugging / standalone tests)
# -------------------------------------------------------------------
if __name__ == "__main__":
    print(repr(Dpiringer()),"\n"*2)
    # static use (without instantiation)
    values = Dpiringer.get_piringer_params("air")
    D = Dpiringer.evaluate("PET",100,40)
    print(f"D1 = {D} [m**2/s]")
    # dynamic use (with instantiation)
    Dmodel = Dpiringer("PET",M=100,T=40)
    Dvalue= Dmodel.eval()
    print(f"D2 = {Dvalue} [m**2/s]")
