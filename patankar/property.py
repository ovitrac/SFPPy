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


@version: 1.21
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2024-07-21
@rev: 2025-03-09

"""

import numpy as np

__all__ = ['ActivityCoefficients', 'DFV', 'Diffusivities', 'Dpiringer', 'Dwelle', 'HenryLikeCoefficients', 'MigrationPropertyModel_validator', 'PartitionCoeffcicients', 'PropertyModelSelector', 'gFHP', 'kFHP', 'migrationProperty']

__project__ = "SFPPy"
__author__ = "Olivier Vitrac"
__copyright__ = "Copyright 2022"
__credits__ = ["Olivier Vitrac"]
__license__ = "MIT"
__maintainer__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@agroparistech.fr"
__version__ = "1.30"

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
    _version = 1.30
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

class HenryLikeCoefficients(migrationProperty):
    property = "Henri-like coefficient"
    notation = "k"
    description = "Mathematical model to estimate Henri-like coefficients"
    SIunits = None

class ActivityCoefficients(migrationProperty):
    property = "Activity coefficient"
    notation = "g"
    description = "Mathematical model to estimate activity coefficients"
    SIunits = None

class PartitionCoeffcicients(migrationProperty):
    property = "Partition Coefficient"
    notation = "K"
    description = "Mathematical model to estimate partition coefficients"
    SIunits = None


# %% Simplified Flory-Huggins model of Henry-like coefficients built on gFHP
class kFHP(HenryLikeCoefficients):
    """
        Simplified model to estimate Henry-like coefficients based on gFHP class
            ki,k = Vi * Pi gik(P'i,P'k,Vi,Vk,crystallinity,porosity)

            i: solute
            k: P or F
            P'i and P'k: Polarity index (e.g.: migrant("solute").polarityindex)
            Vi, Vk: molar volumes (e.g. migrant("solute").molarvolumeMiller)
            ispolymer: True for polymers
            alpha: scaling constant for chiik (default=0.14=1/migrant("water").polarityindex)
            lngmin: minimum value (default=0)
            Psat: vapor saturation pressure
            cristallinity: crystallinity of the solid phase
            porosity: porosity of the effective solid.medium

        Only a static evaluate is proposed.

        Note use: scaling = False to get activity coefficients instead of Henry-like ones

    """
    name = "FHP"
    description = "Flory-Huggins model of Henry-likecoefficients from P' and V at infinite dilution in k"
    model = "semi-empirical"
    theory = "Flory-Huggins"
    parameters = {  "Pi": {"description": "polarity index of solute i","units": "-"},
                    "Pk": {"description": "polarity index of continuous phase k","units": "-"},
                    "Vi": {"description": "molar volume","units": "g/cm**3"},
                    "Vk": {"description": "molar volume of k","units": "g/cm**3"},
                  "Psat": {"description": "vapor saturation pressure of i","units": "Pa"},
         "crystallinity": {"description": "crystallinity of the solid phase",'units':"-"},
              "porosity": {"description": "porosity of the effective solid/medium",'units':"-"},
                }
    _available_to_import = True # this model can be imported

    @classmethod
    def evaluate(cls, Pi=1.41, Pk=3.97, Vi=124.1, Vk=30.9, ispolymer = False, alpha=0.14,lngmin=0.0,Psat=1.0,scaling=True,porosity=0,crystallinity=0):
        """evaluate gFHP model(Pi,Pk,Vi,Vk,ispolymer)"""
        scalesolidamorphous = (1-porosity)*(1-crystallinity)
        scalesolidamorphous = 1 if scalesolidamorphous==0 else scalesolidamorphous # pure air
        if scaling: # (default behavior)
            return gFHP.evaluate(Pi=Pi, Pk=Pk, Vi=Vi, Vk=Vk, ispolymer=ispolymer,
                        alpha=alpha,lngmin=lngmin,
                        gscale=Vi*1e-3*Psat/scalesolidamorphous # Vi is converted [cm**3/g] --> [m**3/kg]
                        )
        else:
            return gFHP.evaluate(Pi=Pi, Pk=Pk, Vi=Vi, Vk=Vk, ispolymer=ispolymer,
                        alpha=alpha,lngmin=lngmin,
                        gscale=1/scalesolidamorphous
                        )

# %% Simplified Flory-Huggins model of activity coefficients
class gFHP(ActivityCoefficients):
    """
        Simplified model to estimate activity coefficients gik from P'i, P'k, Vi, Vk
            i: solute
            k: P or F
            P'i and P'k: Polarity index (e.g.: migrant("solute").polarityindex)
            Vi, Vk: molar volumes (e.g. migrant("solute").molarvolumeMiller)
            ispolymer: True for polymers
            alpha: scaling constant for chiik (default=0.14=1/migrant("water").polarityindex)
            lngmin: minimum value (default=0)
            gscale: activity coefficient (default=1.0)

            Use gscale to enforce gik<=1

        Only a static evaluate is proposed.

    """
    name = "FHP"
    description = "Flory-Huggins model of activity coefficients from P' and V at infinite dilution in k"
    model = "semi-empirical"
    theory = "Flory-Huggins"
    parameters = {"Pi": {"description": "polarity index of solute i","units": "-"},
                  "Pk": {"description": "polarity index of continuous phase k","units": "-"},
                  "Vi": {"description": "molar volume of i","units": "-"},
                  "Vk": {"description": "molar volume of k","units": "-"}
                }
    _available_to_import = True # this model can be imported

    @classmethod
    def evaluate(cls, Pi=1.41, Pk=3.97, Vi=124.1, Vk=30.9, ispolymer = False,
                 alpha=0.14,lngmin=0.0,gscale=1.0):
        """evaluate gFHP model(Pi,Pk,Vi,Vk,ispolymer)"""
        if ispolymer:
            rik = 0.0
            nik = 0.0
            chimin = 0.25
        else:
            rik = Vi/Vk
            nik = (rik - 5)/5
            chimin = 0
        if Pi is None or Pk is None:
            raise RuntimeError("✋🏻🛑⛔️ At least, one of the elements (migrant/medium/polymer) lacks ⌬ chemical information.")
        chiik = np.maximum(chimin,alpha * (Pi - Pk)**2)
        lngik = np.maximum(lngmin,chiik + 1 - (rik - nik))
        return gscale * np.exp(lngik)



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
    parameters = {"polymer":{"polymer": "polymer code/name", "units":"N/A"},
                  "M": {"description": "molecular mass","units": "g/mol"},
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
    # duplicate an entry for wPET from rPET
    piringer_data["wPET"] = piringer_data["rPET"]
    piringer_data["wPET"]["className"] = "wPET"


    # Dpiringer constructor
    def __init__(self, polymer="LDPE", M=100, T=40):
        """
        Instantiate a Dpiringer object for a specific polymer key
        (e.g. 'LDPE', 'PET',...). The corresponding App and tau
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
        return np.exp(exponent)


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
        D = np.exp(exponent)
        return D

# %% DFV model
class DFV(Diffusivities):
    """
        Diffusivity predicted Hole free-volume model from this reference.
        This model covers well plasticizing effects and is applicable for substances built
        on a repeated pattern connecting linearly. Anchor effects are also included.

        Current implementation covers only toluene as surrogate for recycled materials.


        REFERENCE
        Zhu Y., Welle, F. and Vitrac O. A blob model to parameterize polymer hole free volumes and solute diffusion",
        *Soft Matter* **2019**, 15(42), 8912-8932. DOI: https://doi.org/10.1039/C9SM01556F

        ABSTRACT
        Solute diffusion in solid polymers has tremendous applications in packaging,
        reservoir, and biomedical technologies but remains poorly understood. Diffusion
        of non-entangled linear solutes with chemically identical patterns (blobs) deviates
        dramatically in polymers in the solid-state (αlin > 1, Macromolecules 2013, 46, 874)
        from their behaviors in the molten state (αlin = 1, Macromolecules, 2007, 40, 3970).
        This work uses the scale invariance of the diffusivities, D, of linear probes
        D(N·M_blob + M_anchor,T,Tg) = N^(-αlin(T,Tg)) * D(M_blob + M_anchor,T,Tg) comprising
        N identical blobs of mass M_blob and possibly one different terminal pattern (anchor of
        mass M_anchor) to evaluate the amounts of hole-free volume in seven polymers (aliphatic,
        semi-aromatic and aromatic) over a broad range of temperatures (−70 K ≤ T − Tg ≤ 160 K).
        The new parameterization of the concept of hole-free volumes opens the application of
        the free-volume theory (FVT) developed by Vrentas and Duda to practically any polymer,
        regardless of the availability of free-volume parameters. The quality of the estimations
        was tested with various probes including n-alkanes, 1-alcohols, n-alkyl acetates, and
        n-alkylbenzene. The effects of enthalpic and entropic effects of the blobs and the anchor
        were analyzed and quantified. Blind validation of the reformulated FVT was tested
        successfully by predicting from first principles the diffusivities of water and toluene
        in amorphous polyethylene terephthalate from 4 °C to 180 °C and in various other polymers.
        The new blob model would open the rational design of additives with controlled diffusivities
        in thermoplastics.
    """

    name = "FV"
    description = "Hole Free Volume model - current implementation is limited to toluene"
    model = "theory"
    theory = ["free-volume","scaling"]
    parameters = {"polymer":{"polymer": "polymer code/name", "units":"N/A"},
                  "T": {"description": "temperature","units": "degC"},
                  "Tg": {"description": "glass transition temperature","units": "degC"}
                }
    _available_to_import = True # this model can be directly imported

    # Constants
    R = 8.31
    T0K = 273.15 # K
    deltaT = 2   # (K) sharpness of the transition at Tg
    betalin = 1  # Rouse scaling

    # Polymer data (Tg in K) stored in a dictionary.
    _data = {
        'LDPE': {'Tg': 148.15, 'D0': 1.87e-08, 'xi': 0.615,  'ref': 3, 'Ka': 144, 'Kb': 40, 'E': 0, 'r': 0.5},
        'PMMA': {'Tg': 381.15, 'D0': 1.87e-08, 'xi': 0.56,   'ref': 2, 'Ka': 252, 'Kb': 65, 'E': 0, 'r': 0.5},
        'PS':   {'Tg': 373.15, 'D0': 4.8e-08,  'xi': 0.584,  'ref': 2, 'Ka': 144, 'Kb': 40, 'E': 0, 'r': 0.5},
        'PVAc': {'Tg': 305.15, 'D0': 1.87e-08, 'xi': 0.86,   'ref': 4, 'Ka': 142, 'Kb': 40, 'E': 0, 'r': 0.5},
        'gPET': {'Tg': 349.15, 'D0': 1.0205e-08, 'xi': 0.6761, 'ref': 5, 'Ka': 252, 'Kb': 65, 'E': 0, 'r': 0.6153},
        'wPET': {'Tg': 316.15, 'D0': 1.02046e-08, 'xi': 0.6761, 'ref': 5, 'Ka': 252, 'Kb': 65, 'E': 0, 'r': 0.277734},
    }

    # Reference data used to parameterize the polymer (matching ref)
    _references = [
        'Vrentas and Vrentas, 1994',
        'Zielinski and Duda, 1992',
        'Lutzow et al., 1999',
        'Hong, 1995',
        # toluene in PET
        'Welle,2008',
        'Pennarun et al., 2004',
        'Welle,2013',
        'our work (permeation)',
        'our work (sorption)',
    ]


    def __init__(self, polymer="gPET", Tg=76.0, T=40.0):
        """
        Instantiate a DFV object for a specific polymer key
        (e.g. 'LDPE', 'PMMA', or 'PET'). The corresponding D0, xi, Ka, Kb
        are looked up and stored as instance attributes.
        """
        polymer_str = polymer.strip()
        if polymer_str not in self._data:
            raise ValueError(f"No exact match for polymer key: {polymer_str!r}")
        self.polymer = polymer_str
        self.solute = "toluene"
        self.Tg = Tg + self.T0K
        self.Ka = self._lookup("Ka")
        self.Kb = self._lookup("Kb")
        self.D0 = self._lookup("D0")
        self.r = self._lookup("r")
        self.E = self._lookup("E")
        self.xi = self._lookup("xi")

    def _lookup(self,prop):
        """Helper function to lookup a property value from the data dictionary"""
        if prop not in self._data[self.polymer]:
            raise ValueError(f"The property {prop} does not exist for {self.polymer}")
        return self._data[self.polymer][prop]

    def alpha(self,T):
        """alpha for T >= Tg"""
        TK = self.T0K+T
        return 1 + self.Ka / (TK - self.Tg + self.Kb)

    def alphag(self,T):
        """alpha for T < Tg"""
        TK = self.T0K+T
        return 1 + self.Ka / (self.r * (TK - self.Tg) + self.Kb)

    def H(self,T):
        """Heaviside-like function using tanh"""
        TK = self.T0K+T
        return 0.5 * (1 + np.tanh(4 / self.deltaT * (TK - self.Tg)))

    def alphaT(self,T):
        """Composite alpha function that smoothly transitions between alpha and alphag"""
        H = self.H(T)
        return (1-H) * self.alphag(T) + H * self.alpha(T)

    def Plike(self,T):
        """Plike function see publication"""
        return (self.alphaT(T) + self.betalin) / 0.24

    def eval(self,T,**extra):
        """Compute FV D for this polymer"""
        TK = self.T0K+T
        if TK-self.Tg < -self.Kb/self.r + self.deltaT:
            return None # temperature too low for theory at glassy state
        else:
            return self.D0 * np.exp(-self.E / (self.R * TK)) * np.exp(-self.xi * self.Plike(T))

    @classmethod
    def evaluate(cls,polymer="gPET", Tg=76.0, T=40.0, **extra):
        """
        Evaluate D (DFV) for toluene in polymer at T in function of its Tg

        Parameters
        ----------
        polymer : str
            Polymer name (e.g. 'LLDPE', 'LDPE', 'rPET', etc.) as listed in the original data structure.
        T : float
            Temperature in °C, default = 40.
        Tg : float
            Glass Transiton Temperature in °C, default = Tg (PET value).


        Returns
        -------
        float
            The estimated diffusion coefficient in m^2/s of toluene.
        """
        FV = DFV(polymer=polymer,Tg=Tg,T=T)
        return FV.eval(T,**extra)


# %% Welle model
class Dwelle(Diffusivities):
    """
        Diffusivities predicted with the Welle model

        References:

            Ewender J, Welle F. A new method for the prediction of diffusion coefficients in poly(ethylene
            terephthalate)—Validation data. Packag Technol Sci. 2022; 35(5): 405-413.
            https://doi.org:10.1002/pts.2638

            Welle, F. (2021). Diffusion Coefficients and Activation Energies of Diffusion of Organic Molecules
            in Polystyrene below and above Glass Transition Temperature. Polymers, 13(8), 1317.
            https://doi.org/10.3390/polym13081317

    """

    name = "Welle"
    description = "Welle diffusivity model"
    model = "empirical"
    theory = "scaling"
    parameters = {"polymer":{"polymer": "polymer code/name", "units":"N/A"},
                  "T": {"description": "temperature","units": "degC"},
                  "Tg": {"description": "glass transition temperature","units": "degC"},
                  "Vvdw": {"description": "molecular volume 3D","units": "Å³"}
                }
    _available_to_import = True # this model can be directly imported


    # Welle values (the primary key matches the one used in layer)
    welle_data = {
        # a in 1/K, b in cm2/s, c in A3, d in 1/K
        "gPET": {"a": 1.93e-3, "b": 2.27e-6, "c": 11.1, "d":1.50e-4},
        "PS": {"a": 2.59e-3, "b": 7.38e-9, "c": 55.71, "d":2.73e-5},
        "rPS": {"a": 2.44e-3, "b": 6.46e-8, "c": 25.51, "d":7.55e-5}, # rubber PS
        "HIPS": {"a": 2.55e-3, "b": 9.21e-9, "c": 73.28, "d": 2.04e-5},
        "rHIPS": {"a": 2.46e-3, "b": 2.07e-7, "c": 45.00, "d": 2.07e-7}, # rubber HIPS
  # add polymers here
        }

    # Constants
    T0K = 273.15 # K

    def __init__(self, polymer="gPET"):
        """
        Instantiate a Dwelle object for a specific polymer key
        (e.g. 'gPET', 'PS', "rPS", "HIPS", or 'rHIPS'). The corresponding a,b,c,d values
        are looked up and stored as instance attributes.
        """
        polymer_str = polymer.strip()
        if polymer_str not in self.welle_data:
            raise ValueError(f"No exact match for polymer key: {polymer_str!r}")
        self.polymer = polymer_str
        self.a = self._lookup("a")
        self.b = self._lookup("b")
        self.c = self._lookup("c")
        self.d = self._lookup("d")

    def _lookup(self,prop):
        """Helper function to lookup a property value from the welle_data dictionary"""
        if prop not in self.welle_data[self.polymer]:
            raise ValueError(f"The property {prop} does not exist for {self.polymer}")
        return self.welle_data[self.polymer][prop]

    def eval(self,Vvdw,T,**extra):
        """Compute D acoording to the Welle model"""
        TK = self.T0K+T
        return  1.0e-4 * self.b * (Vvdw/self.c) ** ((self.a-1/TK)/self.d) # result in m2/s

    @classmethod
    def evaluate(cls,polymer="gPET", Vvdw=100, T=40.0, **extra):
        """
        Evaluate D (Dwelle) for a substance with a molar volume V in polymer at T

        Parameters
        ----------
        polymer : str
            Polymer name (e.g. 'gPET', 'PS', 'HIPS', 'rPS', 'rHIPS' etc.) as listed in the original data structure.
        Vvdw : float
            3D molecular volume, default = 100 (units in A**3).
        T : float
            Temperature in °C, default = 40.

        Returns
        -------
        float
            The estimated diffusion coefficient in m^2/s
        """
        FW = Dwelle(polymer=polymer)
        return FW.eval(Vvdw,T,**extra)

# %% Available models to expose to layer.py and food.py
# List below the importable models (currently only D, k, K models are possible)
# They are imported via:
#    from property import MigrationPropertyModels, MigrationPropertyModel_validator
# Note that the name of the attribute (eg, "Piringer" or "FV") must match class.name (eg, Dpiringer.name, DFV.name)
# A strict validator is proposed as MigrationPropertyModel_validator()

MigrationPropertyModels = {
    "D":{
        "Piringer": Dpiringer,
        "FV": DFV,
        "Welle": Dwelle,
        # add other diffusivity models here
        },
    "k":{
        "FHP": kFHP
        # add other Henry-like models here
        },
    "g":{
        "FHP": gFHP
        # add other activity coefficients models here
        },
    "K":{
        },
    }

# %% Helper functions

# Function helper to get a strict control on property models used by layer.py and food.py
def MigrationPropertyModel_validator(model=None,name=None,notation=None):
    """ Returns True if the proposed model is valid for the requested migraton property """
    rootclass = "migrationProperty"
    expectedpropclass = {"D":"Diffusivities",
                     "k":"HenryLikeCoefficients",
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



# Function helper to select a model based in rules
def PropertyModelSelector(rules, obj, model1=None, params1=None, model2=None, params2=None, flags=None):
    """
    Selects between two models (and their associated parameter dictionaries) based on a set of rules
    evaluated on a provided object (or objects). New features include optional models/params and
    additional operators.

    ---------------------------------------------------------------------------------------------------
    ===== **Important Notice** =====
    ---------------------------------------------------------------------------------------------------
    Several paradigms are available, it is implemented for patankar.loadpubchem.migrant instances as
            migrant.suggest_alt_Dmodel(material,layerindex) by using global rules: Dmodel_extensions

    This low-level function enforces rigourously complex rules to pick the best model based on objexct
    attributes. Objects can be included in a list to test conditions on material (specific layer),
    substance, etc, all together.

    Please refer to examples and current implementations for details.
    ---------------------------------------------------------------------------------------------------

    Parameters:
        rules (dict or list of dict): A rule or a list of rules. A single rule should have the format:
            {
              "operation": "and" or "or" (default "and"),
              "list": [
                  {"attribute": <str>, "op": <str>, "value": <any>, "index": <int> (optional)},
                  ...
              ]
            }
        obj (dict, object, or list/tuple): An object (or sequence of objects) on which the rules are checked.
            For each condition, the attribute is retrieved either from a dict (via key) or from an object
            (via getattr).
        model1 (function or None): The default model function.
        params1 (dict or None): The parameters for model1.
        model2 (function or None): The alternate model function.
        params2 (dict or None): The parameters for model2.
        flags (dict, optional): A dictionary with flags for string comparisons. Defaults to:
            {
              "remove_blanks": True,      # Remove spaces from the string.
              "trim": True,               # Remove leading/trailing whitespace.
              "case_insensitive": True    # Compare in lowercase.
            }
            These flags are applied to both the attribute value (from obj) and the condition value when
            they are strings.

    Returns:
        - If all of model1, params1, model2, and params2 are None: returns a single boolean value
          (the test result).
        - Otherwise: returns a tuple (testresult, selected_model, selected_params) where:
            * testresult (bool): Outcome of evaluating the rules on obj.
            * If testresult is True, selected_model and selected_params are model2 and params2.
            * Otherwise, they are model1 and params1.

    Note: Calling with a tuple of two objects (migrant, medium):
        result, selected_model, selected_params = PropertyModelSelector(
            Dmodel_extensions["DFV"]["rules"],
            (migrant, medium),
            model1, params1,
            model2, params2
        )

    Advanced features:
    1. Pseudo Recursion for List Inputs:
       If both rules and obj are lists (or tuples), then for each rule in rules, the function applies
       the rule to the corresponding object from obj (if there are more rules than objects, the last
       object is reused). The overall test result is True only if all evaluations are True.
    2. Optional Index Field:
       Each condition may include an optional `index` field. If present and if the attribute value
       is a list or tuple, the condition is evaluated on the element at that index.
    3. Optional Models/Parameters:
       If all of model1, params1, model2, and params2 are None, the function returns a single boolean
       result (the outcome of evaluating the rules on obj). Otherwise, if at least one is provided, the
       function returns a tuple: (testresult, selected_model, selected_params). In that case, model1 and
       params1 must be provided.
    2. Additional Operators:


    List of implemented operators
      - **"=" or "=="**
        Tests for equality between the processed attribute value and the condition value.
      - **"in"**
        Checks whether the processed attribute value is a member of the condition value
        (which can be a string, list, or tuple).
      - **"startswith"**
        For string values, verifies if the processed attribute value starts with the processed
        condition value.
      - **"endwith"**
        For string values, verifies if the processed attribute value ends with the processed
        condition value.
      - **">"**
        Performs a numeric greater-than comparison.
      - **"<"**
        Performs a numeric less-than comparison.
      - **"istrue"**
        Tests whether the attribute value is `True` (the condition's value is optional and
        defaults to `True`).
      - **"isfalse"**
        Tests whether the attribute value is `False`.
      - **"all"**
        Applies Python’s built-in `all()` to the attribute value (expects an iterable).
      - **"any"**
        Applies Python’s built-in `any()` to the attribute value (expects an iterable).
      - **"hasattr"**
        Uses Python’s `hasattr()` to check if the attribute value (which may itself be an object)
        has a specified attribute.
        *(Here, the condition value must be provided as the attribute name.)*
      - **"isinstance"**
        Uses `isinstance(attribute_value, condition_value)` to check if the attribute value is an
        instance of the given type (or tuple of types).
      - **"issubclass"**
        Uses `issubclass(attribute_value, condition_value)` to determine if the attribute value
        (expected to be a class) is a subclass of the specified type, with error handling for
        non-class values.
      - **"callable"**
        Checks if the attribute value is callable (i.e. it is a function, method, or any object
        implementing `__call__`).

    Raises:
        ValueError: If not all models/params are None and model1 or params1 is missing, or if an unsupported
                    operator or invalid operation is encountered.

    Example (for the current implementation, refer to Dmodel_extensions in read patankar.Dmodel_extensions)
        Dmodel_extensions = { # toy example, not applicable for production
            "DFV": {
              "description": "hole Free-Volume theory model for toluene in many polymers",
                  "objects": ["migrant","material"],
                    "rules": [
                                {"list": [{"attribute": "InChiKey",
                                           "op": "==",
                                           "value": "YXFVVABEGXRONW-UHFFFAOYSA-N"}] # <--- migrant must be Toluene
                                },
                                {"list": [
                                    {"attribute": "ispolymer",
                                     "op": "==",
                                     "value": True
                                    }, # <--- medium must be a polymer (ispolymer == True)
                                    {"attribute": "layerclass_history",
                                     "index":0,
                                     "op": "in",
                                     "value": ("gPET","LDPE","PP","PS")
                                    } # <---- medium must be one of these polymers
                                        ]
                                }
                            ]
                    }
            }
        from pprint import pp as disp
        from patankar.loadpubchem import migrant
        from patankar.layer import gPET, PS, PP, LDPE, rigidPVC
        from patankar.property import PropertyModelSelector
        m1 = migrant("toluene")
        m2 = migrant("BHT")
        material = gPET()+PS()+PP()+LDPE()+rigidPVC()
        disp(Dmodel_extensions,depth=7,width=60) # show Dmodel_extensinos

        # check FVT
        # Note that the index Dmodel_extensions["DFV"]["rules"][1]["list"][1]["index"] must be assigned
        mig = m1     # migrant
        index = 2    # layer index
        FVTrules = Dmodel_extensions["DFV"]["rules"].copy()
        FVTrules_layer = FVTrules[1]["list"][1]
        FVTrules_layer["index"] = index # layer index
        print(f"FVT({mig.compound},{FVTrules_layer['value'][index]})=",
              PropertyModelSelector(FVTrules,(mig,material))
              )

    """

    # Check if all model/params are None; if so, we will return a single boolean.
    all_models_none = (model1 is None and params1 is None and model2 is None and params2 is None)

    # Pseudo Recursion: if rules and obj are lists/tuples.
    if isinstance(rules, (list, tuple)) and isinstance(obj, (list, tuple)):
        resbyrules = []
        N = len(obj)
        for i, rule_item in enumerate(rules):
            current_obj = obj[i] if i < N else obj[-1]
            result = PropertyModelSelector(rule_item, current_obj, model1, params1, model2, params2, flags=flags)
            # If returning a tuple, extract the boolean result (first element)
            if isinstance(result, tuple):
                result = result[0]
            resbyrules.append(result)
        final_result = all(resbyrules)
        if all_models_none:
            return final_result
        else:
            return (True, model2, params2) if final_result else (False, model1, params1)

    # If not all models are None, ensure model1 and params1 are provided.
    if not all_models_none:
        if model1 is None or params1 is None:
            raise ValueError("model1 and params1 are required if not all models/params are None")
        if model2 is None and params2 is None:
            return (False, model1, params1)
        if model2 is None or params2 is None:
            raise ValueError("model2 and params2 are required if not all models/params are None")
        if  rules is None: # Only return test result: default is False.
            return (False, model1, params1)
    elif rules is None: # Only return test result: default is False.
            return False

    # Set default flags if not provided.
    if flags is None:
        flags = {"remove_blanks": True, "trim": True, "case_insensitive": True}

    # Helper: Process a string with the provided flags.
    def process_str(s, flags):
        if not isinstance(s, str):
            return s
        if flags.get("trim", True):
            s = s.strip()
        if flags.get("remove_blanks", True):
            s = s.replace(" ", "")
        if flags.get("case_insensitive", True):
            s = s.lower()
        return s

    # Evaluate a single condition.
    def evaluate_condition(condition):
        attr = condition.get("attribute")
        operator = condition.get("op")
        cond_value = condition.get("value", None)  # value is optional for istrue/isfalse
        index = condition.get("index", None)  # Optional index for list/tuple attributes

        # Retrieve the attribute value from obj (dict or object).
        if isinstance(obj, dict):
            if attr not in obj:
                return False
            param_value = obj[attr]
        else:
            if not hasattr(obj, attr):
                return False
            param_value = getattr(obj, attr)

        # If an index is specified and param_value is a list/tuple, select that element.
        if index is not None and isinstance(param_value, (list, tuple)):
            try:
                param_value = param_value[index]
            except IndexError:
                return False

        # For string comparisons, process the values if both are strings.
        if isinstance(param_value, str) and isinstance(cond_value, str):
            proc_param = process_str(param_value, flags)
            proc_cond = process_str(cond_value, flags)
        else:
            proc_param = param_value
            proc_cond = cond_value

        # Standard operators.
        if operator in ("=", "=="):
            return proc_param == proc_cond
        elif operator == "in":
            if isinstance(param_value, str):
                proc_param = process_str(param_value, flags)
                if isinstance(cond_value, str):
                    proc_cond = process_str(cond_value, flags)
                    return proc_param in proc_cond
                elif isinstance(cond_value, (list, tuple)):
                    proc_list = [process_str(item, flags) if isinstance(item, str) else item
                                 for item in cond_value]
                    return proc_param in proc_list
                else:
                    return proc_param in cond_value
            else:
                return proc_param in cond_value
        elif operator == "startswith":
            if not isinstance(param_value, str):
                return False
            return process_str(param_value, flags).startswith(process_str(cond_value, flags))
        elif operator == "endwith":
            if not isinstance(param_value, str):
                return False
            return process_str(param_value, flags).endswith(process_str(cond_value, flags))
        # New operators:
        elif operator == "istrue":
            # value field is optional; equivalent to checking equality with True.
            return proc_param == True
        elif operator == "isfalse":
            return proc_param == False
        elif operator == "all":
            try:
                return all(param_value)
            except Exception:
                return False
        elif operator == "any":
            try:
                return any(param_value)
            except Exception:
                return False
        elif operator == "hasattr":
            # cond_value must be provided as the attribute name to check.
            return hasattr(param_value, cond_value)
        elif operator == "isinstance":
            return isinstance(param_value, cond_value)
        elif operator == "issubclass":
            try:
                return issubclass(param_value, cond_value)
            except Exception:
                return False
        elif operator == "callable":
            return callable(param_value)
        elif operator == ">":
            try:
                return param_value > cond_value
            except Exception:
                return False
        elif operator == "<":
            try:
                return param_value < cond_value
            except Exception:
                return False
        else:
            raise ValueError(f"Unsupported operator: {operator}")

    # Evaluate all conditions from the rule.
    conditions = rules.get("list", [])
    results = [evaluate_condition(cond) for cond in conditions]
    op_str = rules.get("operation", "and").lower()
    if op_str == "and":
        final_result = all(results)
    elif op_str == "or":
        final_result = any(results)
    else:
        raise ValueError("Invalid rules operation. Must be 'and' or 'or'")

    if all_models_none:
        return final_result
    else:
        return (True, model2, params2) if final_result else (False, model1, params1)




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

    # Define dummy functions for testing.
    def dummy1():
        return "dummy1"

    def dummy2():
        return "dummy2"

    # Test with obj as a dictionary.
    test_obj_dict = {"a": " test ", "num": 15}
    rules = {"operation": "and", "list": [
        {"attribute": "a", "op": "==", "value": "Test"},
        {"attribute": "num", "op": ">", "value": 10}
    ]}
    res, m, p = PropertyModelSelector(rules, test_obj_dict, dummy1, {"a": 1}, dummy2, {"a": "alternate"})
    assert res is True and m == dummy2 and p == {"a": "alternate"}, "Expected dummy2 with dict-based obj"

    # Test with obj as a generic object.
    class TestObj:
        def __init__(self, a, num):
            self.a = a
            self.num = num

    test_obj_obj = TestObj(" test ", 15)
    res, m, p = PropertyModelSelector(rules, test_obj_obj, dummy1, {"a": 1}, dummy2, {"a": "alternate"})
    assert res is True and m == dummy2 and p == {"a": "alternate"}, "Expected dummy2 with object-based obj"

    # Test: Attribute missing in object.
    test_obj_obj2 = TestObj(" test ", 5)
    rules_missing = {"operation": "and", "list": [{"attribute": "b", "op": "==", "value": "something"}]}
    res, m, p = PropertyModelSelector(rules_missing, test_obj_obj2, dummy1, {"a": 1}, dummy2, {"a": "alternate"})
    assert res is False and m == dummy1, "Expected default when attribute missing in object"

    print("All tests passed.")


    # other test
    from patankar.loadpubchem import migrant
    m = migrant("toluene")
    istoluene = {"list": [{"attribute": "InChiKey", "op": "==", "value": "YXFVVABEGXRONW-UHFFFAOYSA-N"}]}
    res = PropertyModelSelector(istoluene,m)
