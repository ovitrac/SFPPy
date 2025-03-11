#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""
===============================================================================
SFPPy Module: LoadPubChem
===============================================================================
Retrieves and caches molecular properties from PubChem for use in migration modeling.
Stores both detailed (`full.json`) and lightweight (`simple.json`) records.

**Main Components:**
- **`CompoundIndex`** (Manages the local PubChem database)
    - Caches compound information to avoid repeated queries
    - Uses a local JSON index for fast lookups by name, CAS, or structure
    - Automatically fetches missing compounds from PubChem
- **`migrant`** (Links chemical compounds to migration properties)
    - Retrieves and stores diffusion (`D`) and partition (`K`) coefficients from `property.py`
    - Supports surrogate or theoretical migrants if missing from PubChem
    - Used directly in `layer.py` to define migrating substances in packaging materials

**Integration with SFPPy Modules:**
- Used in `layer.py` to define a migrating chemical with associated properties.
- Calls `property.py` when PubChem does not provide required mass transfer parameters.

Example:
```python
from loadpubchem import migrant
m = migrant(name="anisole")
D_value = m.D.evaluate(polymer="LDPE", T=60)
```


===============================================================================
Details
===============================================================================
This module offers a simple mechanism to retrieve molecular properties from the US-NIH Pubchem website.
The data are cached locally for efficiency and can be used off-line.
The high-level class migrant connects the attributes of the chemical with predictions available in the
module `property.py` including diffusivities and partitioning.


Overview of Compound Index
---------------------------
This module provides a `CompoundIndex` class that allows you to locally cache chemical data retrieved from PubChem, saving each record in two JSON-based formats:

1. **Full JSON** (cidXXXX.full.json) – A comprehensive snapshot of all PubChem properties for a given compound (e.g., synonyms, IUPAC name, molecular weight, 3D data, etc.).
2. **Simple JSON** (cidXXXX.simple.json) – A minimal “lightweight” record containing only the most essential fields, such as:
   - CID, name, synonyms, CAS numbers, molecular weight, formula, SMILES, InChI, InChIKey, and date.

Each record is indexed by any recognized synonym (including CAS, IUPAC names, common names, etc.), so that subsequent searches by any of these identifiers will pull up the same compound. A local synonyms→[CIDs] mapping is stored in a single JSON index file (default: pubchem_index.json).

Key Features
------------
- **Local Cache**: Creates and maintains a folder (default: `cache.PubChem`) where each compound is stored in two forms:
  1. **`cidXXXX.full.json`** – Contains all available properties from PubChem.
  2. **`cidXXXX.simple.json`** – A lighter record used for quick lookups and indexing.

- **Synonym Indexing**: All synonyms (including IUPAC name, title, and anything else treated as a synonym) are captured and mapped to the compound’s CID in a local dictionary, serialized as `pubchem_index.json`. As soon as a new compound is retrieved from PubChem, these synonyms are added to the index so that future searches for any of those synonyms will immediately return the correct record from the local cache.

- **Refreshable Index**: If the index does not exist or is invalid, the module scans all `*.full.json` files in the cache folder, regenerates each `*.simple.json` (if missing or outdated), and rebuilds the synonyms index.

- **PubChem Queries**: When a requested compound is *not* in the local index, the code automatically queries PubChem via your private `pubchempy` library (`get_compounds(...)`). The first match is saved locally, indexed, and returned. If no match is found, an empty result is returned.

- **Flexible Searching**: You can search by compound name, CAS, SMILES, or any other string recognized by PubChem’s "name" lookup. The local index supports direct substring equality (lowercased). You can easily extend or adapt the code for fuzzy matches or partial synonyms.

Usage Example
-------------
Below is a minimal example of how you might use the `CompoundIndex` class once the module is imported:

.. code-block:: python

    from compound_cache import CompoundIndex

    # Instantiate and automatically load or build the index
    db = CompoundIndex()

    # Search for anisole in "simple" mode
    result_simple = db.find("anisole", output_format="simple")
    print("Simple record:\n", result_simple)

    # Retrieve the full record for anisole
    result_full = db.find("anisole", output_format="full")
    print("Full record:\n", result_full)

When you search for a compound, the class:
1. Checks if that query (in lowercase) is in the local synonyms index.
2. If it is found, loads either the `*.simple.json` or `*.full.json` (depending on output_format) and returns one or more results in a pandas DataFrame.
3. If not found, queries PubChem once, saves the new record locally in both file formats, adds synonyms to the index, and returns a one-row DataFrame.

Class Summary
-------------
**`CompoundIndex`**
- **`__init__(cache_dir='cache.PubChem', index_file='pubchem_index.json')`**
  Prepares the local cache folder, loads or refreshes the synonyms index.

- **`refresh_index()`**
  Rebuilds the index by scanning all `*.full.json` files in the local cache. Re-creates any missing `*.simple.json` from the full data.

- **`find(query, output_format='simple')`**
  Main user-facing method. Returns a pandas DataFrame of matching records (possibly multiple rows if synonyms map to multiple CIDs). If no local match is found, it queries PubChem, stores the record, and updates the index.

- **Internal Helper Methods**
  - **`_extract_all_pubchem_properties(compound_obj)`**: Extracts every property from the `pubchempy.Compound` object, calling each property accessor (cid, synonyms, iupac_name, etc.).
  - **`_generate_simple_dict(full_data, synonyms_set=None)`**: Produces the minimal “light” dictionary saved to `cidXXXX.simple.json`.
  - **`_gather_synonyms(full_data)`**: Merges synonyms (and other text fields, like `iupac_name`/`title`) into a unified set of strings.
  - **`_add_synonym_to_index(synonym, cid)`**: Inserts or updates the synonyms→[CIDs] mapping.

Dependencies
------------
- **`pandas`**: For returning results as DataFrame objects.
- **`json`** & **`os`** & **`glob`** & **`datetime`**: For file I/O, directory handling, indexing, caching.
- **`re`**: For simple CAS pattern matching (e.g., `^\d{1,7}-\d{2}-\d$`).
- **`private.pubchempy`**: A local (private) version of the PubChem Python client, providing the `get_compounds(...)` method and `Compound` property accessors.

Note
-----
- Large-scale usage: For large compound sets, consider optimizing the index or storing data in a more robust database.
- The synonyms approach: Default matching is **exact** (lowercased). Fuzzy or partial matches require custom logic.


@version: 1.30
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2024-03-10
@rev: 2025-03-06

Version History
---------------
- 1.0: Initial version, supporting local caching, synonyms index, and direct PubChem lookup.
- 1.2: Production
- 1.21: PubChem cap rate enforced (urgent request)

"""


import os
import subprocess
import requests
import json
import re
import glob
import pandas as pd
import numpy as np
from datetime import datetime
import time

try:
    from PIL import Image
    PIL_AVAILABLE = True
except ImportError:
    PIL_AVAILABLE = False

# private version of pubchempy
from patankar.private.pubchempy import get_compounds

__all__ = ['CompoundIndex', 'dbdefault', 'get_compounds', 'migrant', 'migrantToxtree', 'polarity_index']

__project__ = "SFPPy"
__author__ = "Olivier Vitrac"
__copyright__ = "Copyright 2022"
__credits__ = ["Olivier Vitrac"]
__license__ = "MIT"
__maintainer__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@agroparistech.fr"
__version__ = "1.30"

# %% Private functions and constants (used by estimators)

# full path of patankar/ used cache.PubChem, cache.Toxtree, private/toxtree/
_PATANKAR_FOLDER = os.path.dirname(__file__)

# Enforcing rate limiting cap: https://www.ncbi.nlm.nih.gov/books/NBK25497/
PubChem_MIN_DELAY = 1 / 3.0  # 1/3 second (333ms)
PubChem_lastQueryTime = 0 # global variable

# returns polarity index from logP and V
def polarity_index(logP=None, V=None, name=None,
                   Vw=19.588376948550433,  # migrant("water").molarvolumeMiller
                   Vo=150.26143432234372,  # migrant("octanol").molarvolumeMiller
                   A=0.18161296829146106,
                   B=-3.412678660396018,
                   C=14.813767205916765):
    """
    Computes the polarity index (P') from a given logP value and molar volume V.
    This is done using a quadratic model fitted to experimental data:

        E = A * (P')² + B * P' + C
        P' = (-B - sqrt(B² - 4A(C - E))) / (2A)

        where:
        - E = logP * ln(10) - S = Xw - Xo
        - S = entropy contribution = - (V/Vw - V/Vo)

    Parameters
    ----------
    logP : float, list, or np.ndarray
        The logP value(s) for which to compute the polarity index P'.
    V : float, list, or np.ndarray
        The molar volume(s) corresponding to logP. Must be either:
        - The same size as `logP`, or
        - A single scalar value that will be applied to all logP values.
    name : str, optional
        A solvent name (instead of providing logP and V). If given, logP and V
        will be fetched from the `migrant` database.
    Vw : float, optional
        Molar volume of water (default: 19.59).
    Vo : float, optional
        Molar volume of octanol (default: 150.26).
    A, B, C : float, optional
        Coefficients for the quadratic equation.

    Returns
    -------
    float or np.ndarray
        The calculated polarity index P'. If logP is out of the valid range:
        - Returns **10.2** for very polar solvents (beyond water).
        - Returns **0** for extremely hydrophobic solvents (beyond n-Hexane).

    Raises
    ------
    ValueError
        If both `logP` and `V` are not provided, or if their lengths do not match.

    Example Usage
    -------------
    >>> logP = migrant("anisole").logP
    >>> V = migrant("anisole").molarvolumeMiller
    >>> polarity_index(logP, V)
    8.34  # Example output

    >>> polarity_index(logP=[-1.0, 0.5, 2.0], V=50)
    array([9.2, 4.5, 1.8])  # Example outputs
    """

    # Define valid logP range based on quadratic model limits
    Emin = C - B**2 / (4*A)  # ≈ -2.78 (theoretical minimum lnKow=E)
    Emax = C                 # ≈ 14.81 (theoretical maximum logP)
    Pmax = 10.2  # Saturation value for highly polar solvents

    # Fetch logP and V if `name` is given
    if logP is None or V is None:
        if name is None:
            raise ValueError("Provide either (logP, V) pair or a valid solvent name.")
        from patankar.loadpubchem import migrant
        tmp = migrant(name)
        logP, V = tmp.logP, tmp.molarvolumeMiller

    # Convert inputs to NumPy arrays for consistency
    logP = np.asarray(logP, dtype=np.float64)
    if np.isscalar(V):
        V = np.full_like(logP, V, dtype=np.float64)  # Broadcast scalar V
    else:
        V = np.asarray(V, dtype=np.float64)

    # Ensure logP and V have matching sizes
    if logP.shape != V.shape:
        raise ValueError("logP and V must have the same shape or V must be a scalar.")

    def compute_P(logP_value, V_value):
        """Computes P' for a single logP and V value after input validation."""
        S = - (1/Vw - 1/Vo) * V_value
        E = logP_value * 2.302585092994046 - S  # Convert logP to natural log (ln)

        # Handle extreme values
        if E < Emin:
            return Pmax  # Very polar solvents
        if E > Emax:
            return 0.0  # Extremely hydrophobic solvents

        # Solve quadratic equation
        discriminant = B**2 - 4*A*(C - E)
        sqrt_discriminant = np.sqrt(discriminant)
        P2root = (-B - sqrt_discriminant) / (2*A)  # Always select P2

        return P2root if P2root <= Pmax else Pmax

    # Vectorized computation for arrays
    return np.vectorize(compute_P)(logP, V)



# %% Core class (low-level)
class CompoundIndex:
    """
    Class to query chemical compounds by any synonym (name, CAS, etc.)
    using a local PubChem cache, refreshing or populating it automatically
    from actual PubChem queries if needed.
    """

    def __init__(self, cache_dir="cache.PubChem", index_file="pubchem_index.json"):
        """
        Constructor: ensures cache directory and index file exist/are loaded.

        :param cache_dir: path to local cache of *.json files
        :param index_file: local JSON file holding synonyms → [cids] index
        """
        self.cache_dir = os.path.join(_PATANKAR_FOLDER,cache_dir)
        os.makedirs(self.cache_dir, exist_ok=True)

        self.index_file = os.path.join(self.cache_dir, index_file)
        # Regex to identify CAS-like strings, e.g. "1234-56-7"
        self._cas_regex = re.compile(r'^\d{1,7}-\d{2}-\d$')

        # Attempt to load existing index; if missing or invalid, rebuild
        if not os.path.isfile(self.index_file):
            self.refresh_index()
        else:
            with open(self.index_file, "r", encoding="utf-8") as f:
                try:
                    # for debugging, we split reading and parsing
                    #rawjson = f.read()
                    #self.index = json.loads(rawjson)
                    self.index = json.load(f)
                except json.JSONDecodeError:
                    print("LOADPUBCHEM: JSON ERROR in {self.index_file}, the current index is discarded.")
                    self.index = {}
            if not isinstance(self.index, dict) or not self.index:
                self.refresh_index()

    def refresh_index(self):
        """
        Rebuild the synonyms→[cids] index by scanning *.full.json files
        in the cache directory, and regenerating each *.simple.json if needed.
        """
        self.index = {}
        full_files = glob.glob(os.path.join(self.cache_dir, "cid*.full.json"))

        for full_path in full_files:
            filename = os.path.basename(full_path)  # e.g. "cid12345.full.json"
            cid_str = filename.replace("cid", "").replace(".full.json", "")
            try:
                cid = int(cid_str)
            except ValueError:
                continue  # skip any weirdly named files

            # Load full data
            with open(full_path, "r", encoding="utf-8") as f:
                try:
                    full_data = json.load(f)
                except:
                    continue

            # Gather synonyms from the "full" data
            synonyms_set = self._gather_synonyms(full_data)

            # Possibly regenerate the *.simple.json
            simple_dict = self._generate_simple_dict(full_data, synonyms_set)
            simple_path = os.path.join(self.cache_dir, f"cid{cid}.simple.json")
            with open(simple_path, "w", encoding="utf-8") as fw:
                json.dump(simple_dict, fw, indent=2)

            # Add synonyms to the index
            for syn in simple_dict.get("synonyms", []):
                self._add_synonym_to_index(syn, cid)

        # Save updated index
        with open(self.index_file, "w", encoding="utf-8") as f:
            json.dump(self.index, f, indent=2)

    def _add_synonym_to_index(self, synonym, cid):
        """
        Helper to map a single synonym→cid in self.index.
        """
        syn_lower = synonym.strip().lower()
        if syn_lower not in self.index:
            self.index[syn_lower] = []
        if cid not in self.index[syn_lower]:
            self.index[syn_lower].append(cid)

    def _gather_synonyms(self, full_data):
        """
        Gathers synonyms from the loaded full-data dictionary.
        We expect 'synonyms' to be a list, plus possible extra fields.
        Merge them into a single set for deduplication.
        """
        synonyms_set = set()

        # If your full_data includes a 'synonyms' list
        syn_list = full_data.get("synonyms", [])
        if syn_list:
            synonyms_set.update(syn_list)

        # Also merge other textual fields that are effectively synonyms
        iupac_name = full_data.get("iupac_name")
        if iupac_name:
            synonyms_set.add(iupac_name)

        title = full_data.get("title")
        if title:
            synonyms_set.add(title)

        # You can add more fields if you treat them as synonyms or common names
        return synonyms_set

    def _generate_simple_dict(self, full_data, synonyms_set=None):
        """
        Builds a small "light" dictionary for quick searching:
            CID, name, synonyms, CAS, M, formula, SMILES, InChi, InChiKey, logP, date.
        """
        if synonyms_set is None:
            synonyms_set = set()

        cid = full_data.get("cid", None)
        synonyms_list = sorted(s.strip() for s in synonyms_set if s.strip())

        # Identify CAS numbers within synonyms
        cas_list = []
        for syn in synonyms_list:
            if self._cas_regex.match(syn):
                cas_list.append(syn)

        # Derive a main 'name'
        name = full_data.get("iupac_name") or (synonyms_list[0] if synonyms_list else "")

        # Some fields might be missing or None
        record = {
            "CID": cid,
            "name": name,
            "synonyms": synonyms_list,
            "CAS": cas_list,
            "M": float(full_data.get("molecular_weight")),
            "formula": full_data.get("molecular_formula"),
            "SMILES": full_data.get("canonical_smiles"),
            "InChi": full_data.get("inchi"),
            "InChiKey": full_data.get("inchikey"),
            "logP": float(full_data.get("xlogp")),
            "date": datetime.now().strftime("%Y-%m-%d"),
        }
        return record

    def _extract_all_pubchem_properties(self, compound_obj):
        """
        Uses your local pubchempy.Compound’s @property accessors to get
        all available fields. This replicates the entire set of property
        definitions you shared (cid, synonyms, iupac_name, xlogp, etc.),
        then returns them in one dict.

        We'll read each property from the compound_obj and store it.
        If your code snippet has more 3D property calls, just do the same.
        """
        d = {}

        # Basic identifiers
        d["cid"] = compound_obj.cid
        # synonyms is a memoized_property
        # so if we do compound_obj.synonyms, it triggers an extra request for synonyms
        d["synonyms"] = compound_obj.synonyms or []

        # Extract “static” properties
        # Many parse data from compound_obj.record['props'] or similar.
        d["sids"] = compound_obj.sids or []
        d["aids"] = compound_obj.aids or []
        d["elements"] = compound_obj.elements
        d["atoms"] = [self._atom_to_dict(a) for a in compound_obj.atoms]  # or just store them raw
        d["bonds"] = [self._bond_to_dict(b) for b in compound_obj.bonds]
        d["coordinate_type"] = compound_obj.coordinate_type
        d["charge"] = compound_obj.charge
        d["molecular_formula"] = compound_obj.molecular_formula
        d["molecular_weight"] = compound_obj.molecular_weight
        d["canonical_smiles"] = compound_obj.canonical_smiles
        d["isomeric_smiles"] = compound_obj.isomeric_smiles
        d["inchi"] = compound_obj.inchi
        d["inchikey"] = compound_obj.inchikey
        d["iupac_name"] = compound_obj.iupac_name
        d["xlogp"] = compound_obj.xlogp
        d["exact_mass"] = compound_obj.exact_mass
        d["monoisotopic_mass"] = compound_obj.monoisotopic_mass
        d["tpsa"] = compound_obj.tpsa
        d["complexity"] = compound_obj.complexity
        d["h_bond_donor_count"] = compound_obj.h_bond_donor_count
        d["h_bond_acceptor_count"] = compound_obj.h_bond_acceptor_count
        d["rotatable_bond_count"] = compound_obj.rotatable_bond_count
        d["fingerprint"] = compound_obj.fingerprint
        # cactvs_fingerprint might be large but let's store it
        d["cactvs_fingerprint"] = compound_obj.cactvs_fingerprint
        d["heavy_atom_count"] = compound_obj.heavy_atom_count
        d["isotope_atom_count"] = compound_obj.isotope_atom_count
        d["atom_stereo_count"] = compound_obj.atom_stereo_count
        d["defined_atom_stereo_count"] = compound_obj.defined_atom_stereo_count
        d["undefined_atom_stereo_count"] = compound_obj.undefined_atom_stereo_count
        d["bond_stereo_count"] = compound_obj.bond_stereo_count
        d["defined_bond_stereo_count"] = compound_obj.defined_bond_stereo_count
        d["undefined_bond_stereo_count"] = compound_obj.undefined_bond_stereo_count
        d["covalent_unit_count"] = compound_obj.covalent_unit_count

        # 3D data (if present)
        d["volume_3d"] = compound_obj.volume_3d
        d["multipoles_3d"] = compound_obj.multipoles_3d
        d["conformer_rmsd_3d"] = compound_obj.conformer_rmsd_3d
        d["effective_rotor_count_3d"] = compound_obj.effective_rotor_count_3d
        d["pharmacophore_features_3d"] = compound_obj.pharmacophore_features_3d
        d["mmff94_partial_charges_3d"] = compound_obj.mmff94_partial_charges_3d
        d["mmff94_energy_3d"] = compound_obj.mmff94_energy_3d
        d["conformer_id_3d"] = compound_obj.conformer_id_3d
        d["shape_selfoverlap_3d"] = compound_obj.shape_selfoverlap_3d
        d["feature_selfoverlap_3d"] = compound_obj.feature_selfoverlap_3d
        d["shape_fingerprint_3d"] = compound_obj.shape_fingerprint_3d

        return d

    def _atom_to_dict(self, atom_obj):
        """
        Optional: convert a pubchempy.Atom instance to a small dict
        with (aid, element, x, y, z, charge, ...).
        """
        return {
            "aid": atom_obj.aid,
            "element": atom_obj.element,
            "x": atom_obj.x,
            "y": atom_obj.y,
            "z": atom_obj.z,
            "charge": atom_obj.charge,
        }

    def _bond_to_dict(self, bond_obj):
        """
        Optional: convert a pubchempy.Bond instance to a small dict
        with (aid1, aid2, order, etc.).
        """
        return {
            "aid1": bond_obj.aid1,
            "aid2": bond_obj.aid2,
            "order": bond_obj.order,
        }

    def find(self, query, output_format="simple"):
        """
        Main method to find a compound from local index or from PubChem.
        Returns a pd.DataFrame with matching records. If multiple CIDs
        match that synonym, returns multiple rows.

        :param query: string synonym/identifier (name, CAS, SMILES, etc.)
        :param output_format: 'simple' or 'full'
        :return: pd.DataFrame with the results (possibly multiple rows)
        """
        global PubChem_lastQueryTime
        qlower = query.strip().lower()

        if qlower not in self.index:
            # Not found locally => do a PubChem call while respecting cap limit
            # doing more than 3 queries per second will ban you for a day or so
            elapsed = time.time() - PubChem_lastQueryTime # time elapsed since last request
            if elapsed < PubChem_MIN_DELAY:
                wait_time = PubChem_MIN_DELAY - elapsed
                print(f"LOADPUBCHEM: Rate limit reached. Waiting {wait_time:.2f} s...")
                time.sleep(wait_time)
            matches = get_compounds(query, 'name')
            PubChem_lastQueryTime - time.time() # update last request time
            if not matches:
                return pd.DataFrame()  # no hits at all

            best = matches[0]
            # Build the "all-props" dictionary from the pubchempy.Compound
            best_dict = self._extract_all_pubchem_properties(best)

            cid = best_dict.get("cid", None)
            if cid is None:
                return pd.DataFrame()  # some edge case with no cid

            # Save the "full" record
            full_name = f"cid{cid}.full.json"
            full_path = os.path.join(self.cache_dir, full_name)
            with open(full_path, "w", encoding="utf-8") as fw:
                json.dump(best_dict, fw, indent=2)

            # Now prepare the synonyms set from that new record
            synonyms_set = self._gather_synonyms(best_dict)
            # Generate the "simple" record
            simple_dict = self._generate_simple_dict(best_dict, synonyms_set)

            # Save the "simple" record
            simple_name = f"cid{cid}.simple.json"
            simple_path = os.path.join(self.cache_dir, simple_name)
            with open(simple_path, "w", encoding="utf-8") as fw:
                json.dump(simple_dict, fw, indent=2)

            # Update the index with synonyms
            for syn in simple_dict.get("synonyms", []):
                self._add_synonym_to_index(syn, cid)
            # Also index the raw query itself
            self._add_synonym_to_index(query, cid)

            with open(self.index_file, "w", encoding="utf-8") as f:
                json.dump(self.index, f, indent=2)

            # Return a single-row DataFrame
            if output_format == "full":
                return pd.DataFrame([best_dict])
            else:
                return pd.DataFrame([simple_dict])

        else:
            # Found in local index => load data from cache
            cids = self.index[qlower]
            results = []
            for cid in cids:
                if output_format == "full":
                    fpath = os.path.join(self.cache_dir, f"cid{cid}.full.json")
                else:
                    fpath = os.path.join(self.cache_dir, f"cid{cid}.simple.json")

                if not os.path.isfile(fpath):
                    continue  # skip if missing or corrupted
                with open(fpath, "r", encoding="utf-8") as f:
                    data = json.load(f)
                results.append(data)

            if not results:
                return pd.DataFrame()

            return pd.DataFrame(results)

# %% Configuration of migrant class: cache, available Dmodel and kmodel

# Main compound database
dbdefault = CompoundIndex(cache_dir="cache.PubChem", index_file="pubchem_index.json")

# Model extensions to be tested with PropertyModelSelector()
"""
# ==============================================================================================
#  Use this emplate to validate the validation of alternative models by PropertyModelSelector()
#    The syntax is sophisticated and accept multiple paradims and criteria.
#
#    see patankar.property.PropertyModelSelector documentation for more details
# ==============================================================================================

# Evaluate first Dmodel_extensions

# Dependencies
from pprint import pp as disp
from patankar.loadpubchem import migrant
from patankar.layer import gPET, PS, PP, LDPE, rigidPVC
from patankar.property import PropertyModelSelector

# Case study
m1 = migrant("toluene")
m2 = migrant("BHT")
material = gPET()+PS()+PP()+LDPE()+rigidPVC()
disp(Dmodel_extensions,depth=7,width=60) # show Dmodel_extensinos

# we build objects on which rules will be tested
# use material first since all Dmodels check against material, but not all against migrant
# index refers the the layer index (it should be last)
objects = {"material":material,"migrant":m1,"index":2}

# check Dmodels with objects
availableExtensions = list(Dmodel_extensions.keys())
applicableExtensions = [False]*len(availableExtensions)
for imodel,modelCode in enumerate(availableExtensions):
    modelObjects = [objects[o] for o in Dmodel_extensions[modelCode]["objects"]]
    modelRules = Dmodel_extensions[modelCode]["rules"].copy()
    modelRules_layer = modelRules[0]["list"][1] # 0=polymer rules, 1=polymer type
    modelRules_layer["index"] = objects["index"] # layer index
    applicableExtensions[imodel] = PropertyModelSelector(modelRules,modelObjects)
    # nice display
    print(f'{modelCode}('
          f'{objects["material"][objects["index"]].layerclass_history[0]},'
          f'{objects["migrant"].compound}) = ',
          applicableExtensions[imodel]
          )
first_true_index = next((i for i, val in enumerate(applicableExtensions) if val), None)
if first_true_index is not None:
    firstApplicableExtension = availableExtensions[first_true_index]
    print(f"first applicable model: {firstApplicableExtension}")
else:
    print("no alternative model, use the default one")
"""


# Alternative Dmodels
Dmodel_extensions = {
    "DFV": {
      "description": "hole Free-Volume theory model for toluene in many polymers",
          "objects": ["material","migrant"], # it requires material and migrant (material always first)
            "rules": [ # we assume AND between all conditions
                        # Condition on the material (always first)
                        {"list": [
                            # medium must be a polymer (ispolymer == True)
                            {"attribute": "ispolymer",
                             "op": "istrue",
                            },
                            # the polymer with index must be of these types
                            {"attribute": "layerclass_history",
                             "index":0,
                             "op": "in",
                             "value": ("gPET","wPET","PMMA","PS","PVAc","LDPE")
                            }, # next condition
                                ] # close list of rules for rules[0]
                        }, # close rules[0]
                        # Condition on migrant
                        {"list": [
                            # migrant must be Toluene (based on its InChiKey)
                            {"attribute": "InChiKey",
                                   "op": "==",
                                   "value": "YXFVVABEGXRONW-UHFFFAOYSA-N"
                            }
                                ], # next condition
                        }, # next rule
                    ] # close rules
            }, # next model

    "Dwelle": {
    "description": "Frank-Welle diffusivity model based on VdW volumes",
        "objects": ["material"],
          "rules": [ # we assume AND between all conditions
                    # Condition on the material (always first)
                    {"list": [
                        # medium must be a polymer (ispolymer == True)
                        {"attribute": "ispolymer",
                         "op": "istrue",
                        },
                        # the polymer with index must be of these types
                        {"attribute": "layerclass_history",
                         "index":0,
                         "op": "in",
                         "value": ("gPET","PS","rPS","HIPS","rHIPS")
                        }, # next condition
                            ] # close list of rules for rules[0]
                    }, # close rules[0]
                    ] # close rules
            }
    }

# Alternative kmodel
kmodel_extensions = {} # they are not implemented yet

# %% Class migrant (high-level)
# =========================================
#               Migrant class
# =========================================

class migrant:
    """
    A class representing a migrating chemical substance.

    It can be initialized in three main ways:

    1) Case (a) - By a textual name/CAS only (for a real compound search):
       ---------------------------------------------------------
       Example:
           m = migrant(name="anisole", db=my_compound_index)
           # or
           m = migrant(name="anisole", db=my_compound_index, M=None)
       In this mode:
         • A lookup is performed using db.find(name), which may return one or more records.
         • If multiple records match, data from each record is merged:
             - compound  = The text used in the query (e.g. "anisole")
             - name      = Concatenation of all distinct names from the search results
             - CAS       = Concatenation of all CAS numbers from the search results
             - M         = The minimum of all found molecular weights, stored in self.M (a numpy array also keeps the full set)
             - formula   = The first formula
             - logP      = All logP values concatenated into a numpy array (self.logP_array).
                           The main attribute self.logP will be the same array or you may pick a single representative.

    2) Case (b) - By numeric molecular weight(s) alone (generic substance):
       ---------------------------------------------------------
       Example:
           m = migrant(M=200)
           m = migrant(M=[100, 500])  # Possibly a range
       In this mode:
         • No search is performed.
         • name = "generic" (unless you override it).
         • compound = "single molecular weight" if 1 entry in M, or
                      "list of molecular weights ranging from X to Y" if multiple.
         • CAS = None
         • M   = the minimum of all provided M values (also stored in a numpy array)
         • logP = None by default, or can be supplied explicitly as an array

    3) Case (c) - Name + numeric M/logP => Surrogate / hypothetical:
       ---------------------------------------------------------
       Example:
           m = migrant(name="mySurrogate", M=[200, 250], logP=[2.5, 3.0])
         or
           m = migrant(name="surrogate", M=200)
       In this mode:
         • No lookup is performed. This is a “fake” compound not found in PubChem.
         • compound = "single molecular weight" or
                      "list of molecular weights ranging from X to Y" if multiple.
         • name = whatever user provides
         • CAS = None
         • M   = min of the provided M array, stored in a numpy array
         • logP = user-provided array or single float, stored in a numpy array

    Attributes
    ----------
    compound : str
        For case (a) => the search text;
        For case (b,c) => textual description of the numeric M array.
    name : str or list
        For case (a) => aggregated list of all found names (string-joined);
        For case (b) => "generic" or user-supplied name;
        For case (c) => user-supplied name.
    CAS : list or None
        For case (a) => aggregated CAS from search results;
        For case (b,c) => None.
    M : float
        The *minimum* M from either the search results or the user-supplied array.
    M_array : numpy.ndarray
        The full array of all M values found or provided.
    logP : float or numpy.ndarray or None
        For case (a) => an array of all logP from the search results (or None if not found);
        For case (b) => None or user-supplied value/array;
        For case (c) => user-supplied value/array.
    """

    # class attribute, maximum width
    _maxdisplay = 40

    # migrant constructor with its predictor templates
    def __init__(self, name=None,   # substance identified by name
                 M=None, logP=None, # substance identified by M, logP(less reliable)

                 Dmodel = "Piringer", # <--- default D model

                 Dtemplate = {"polymer":"LLDPE",
                              "M":50.0,  # used by Dpiringer (molecular mass in g/mol)
                           "Vvdw":100.0, # used by Dwelle (molecular volume)
                              "T":40.0,  # used by Dpringer, DFV
                             "Tg":76.0,  # used by DFV
                              }, # do not use None

                 kmodel = "FHP",    # <--- default k model

                 ktemplate = {"Pi":1.41,  # P'i (polarity index)
                              "Pk":3.97,  # P'k (polarity index)
                              "Vi":124.1, # molar volume of i
                              "Vk":30.9,  # molar volume of k
                              "ispolymer":True, # True if FH theory is applicable
                              "alpha":0.14, # \alpha \times (P'_i-P'k)^2
                              "lngmin":0.0, # min of log(\gamma_i) -- see theory at inifinite dilution
                              "Psat":1.0,   # partial saturation pressure for i (usually not defined)
                              "crystallinity":0, # k are calculated respectively to the volume fraction
                              "porosity":0       # of amorphous phase (1-crystallinity)(1-porosity)
                              }, # do not use None

                 db=dbdefault, # cache.PubChem database

                 raiseerror=True # raise an error if the susbtance is not found

                 ):
        """
        Create a new migrant instance.

        Parameters
        ----------
        name : str or None
            - A textual name for the substance to be looked up in PubChem (case a),
              or a custom name for a surrogate (case c).
            - If None, and M is given, we treat it as a numeric-only initialization (case b).
        M : float or list/ndarray of float or None
            - For case (a): If provided as None, we do a PubChem search by name.
            - For case (b): The numeric molecular weight(s). No search is performed if name is None.
            - For case (c): Combined name and numeric M => a surrogate with no search.
        logP : float or list/ndarray of float or None
            - For case (a): Typically None. If the PubChem search returns logP, it’s stored automatically.
            - For case (b,c): user can supply. If given, stored in self.logP as a numpy array.
        db : instance of CompoundIndex or similar, optional
            - If you want to perform a PubChem search (case a) automatically, pass an instance.
            - If omitted or None, no search is attempted, even if name is given.
        raiseerror : bool (default=True), optional
            Raise an error if name is not found

        Advanced Parameters
        -------------------
        Property models from MigrationPropertyModels can be directly attached to the substance.
        Based on the current version of migration.py two models are proposed:
            - Set a diffusivity model using
                    - Dmodel="model name"
                      default ="Piringer"
                    - Dtemplate=template dict coding for the key:value parameters
                      (e.g, to bed used Diringer(key1=value1...))
                      note: the template needs to be valid (do not use None)
                      default = {"polymer":None, "M":None, "T":None}
            - Set a Henry-like model using
                    - kmodel="model name"
                      default =None
                    - ktemplate=template dict coding for the key:value parameters
                      default =  {"Pi":1.41, "Pk":3.97, "Vi":124.1, "Vk":30.9, "ispolymer":True, "alpha":0.14, "lngmin":0.0,"Psat":1.0}
            other models could be implemented in the future, read the module property.py for details.

        Example of usage of Dpiringer
            m = migrant(name='limonene')
            # without the helper function
            Dvalue = m.D.evaluate(**dict(m.Dtemplate,polymer="LDPE",T=60))
            # with the helper function
            Dvalue = m.Deval(polymer="LDPE",T=60)

        Raises
        ------
        ValueError if insufficient arguments are provided for any scenario.
        """

        # local import
        # import implicity property migration models (e.g., Dpiringer)
        from patankar.property import MigrationPropertyModels, MigrationPropertyModel_validator

        self.compound = None   # str
        self.name = None       # str or list
        self.cid = None        # int or list
        self.CAS = None        # list or None
        self.M = None          # float
        self.formula = None
        self.smiles = None
        self.M_array = None    # np.ndarray
        self.logP = None       # float / np.ndarray / None

        # special case
        if name==M==None:
            name = 'toluene'

        # Convert M to a numpy array if given
        if M is not None:
            if isinstance(M, (float, int)):
                M_array = np.array([float(M)], dtype=float)
            else:
                # Convert to array
                M_array = np.array(M, dtype=float)
        else:
            M_array = None

        # Similarly, convert logP to array if provided
        if logP is not None:
            if isinstance(logP, (float, int)):
                logP_array = np.array([float(logP)], dtype=float)
            else:
                logP_array = np.array(logP, dtype=float)
        else:
            logP_array = None

        # Case (a): name is provided, M=None => real compound lookup
        if (name is not None) and (M is None):
            if db is None:
                raise ValueError("A db instance is required for searching by name when M is None.")

            df = db.find(name, output_format="simple")
            if df.empty:
                if raiseerror:
                    raise ValueError(f"<{name}> not found")
                print(f"LOADPUBCHEM ERRROR: <{name}> not found - empty object returned")
                self.compound = name
                self.name = [name]
                self.cid = []
                self.CAS = []
                self.M_array = np.array([], dtype=float)
                self.M = None
                self.formula = None
                self.smiles = None
                self.logP = None
            else:
                self.compound = name
                # Initialize dictionary to accumulate properties from each row
                all_data = {
                    "names": [],
                    "cid": [],
                    "cas": [],
                    "inchi": [],
                    "inchikey": [],
                    "m": [],
                    "logp": [],
                    "formula": [],
                    "smiles": [],
                }

                # Helper functions
                def ensure_list(val):
                    if val is None:
                        return []
                    return val if isinstance(val, list) else [val]

                def to_float(val):
                    try:
                        return float(val)
                    except Exception:
                        return np.nan

                for _, row in df.iterrows():
                    # Combine name and synonyms into one set
                    names = ensure_list(row.get("name", []))
                    syns = ensure_list(row.get("synonyms", []))
                    all_data["names"].extend(set(names) | set(syns))

                    # Process CID (append if exists)
                    cid = row.get("CID", None)
                    if cid:
                        all_data["cid"].append(cid)

                    # Process CAS (extend list if exists)
                    cas = ensure_list(row.get("CAS", []))
                    all_data["cas"].extend(cas)

                    # Process InChi and InChiKey (ensure list conversion and extend)
                    inchi = ensure_list(row.get("InChi", []))
                    inchikey = ensure_list(row.get("InChiKey", []))
                    all_data["inchi"].extend(inchi)
                    all_data["inchikey"].extend(inchikey)

                    # Process M: convert to float or use np.nan
                    m_val = row.get("M", None)
                    all_data["m"].append(to_float(m_val) if m_val is not None else np.nan)

                    # Process logP similarly (if not None or empty)
                    logp_val = row.get("logP", None)
                    if logp_val not in (None, ""):
                        all_data["logp"].append(to_float(logp_val))
                    else:
                        all_data["logp"].append(np.nan)

                    # Process formula and SMILES (append even if None to preserve index)
                    all_data["formula"].append(row.get("formula", None))
                    all_data["smiles"].append(row.get("SMILES", None))

                # Convert lists to arrays for numerical properties
                arr_m = np.array(all_data["m"], dtype=float)
                arr_logp = np.array(all_data["logp"], dtype=float)

                # Deduplicate fields where necessary
                unique_names = list(set(all_data["names"]))
                unique_cid = list(set(all_data["cid"]))
                unique_cas = list(set(all_data["cas"]))
                unique_inchi = list(set(all_data["inchi"]))
                unique_inchikey = list(set(all_data["inchikey"]))

                # Store deduplicated and processed values in the object
                self.name = unique_names
                self.cid = unique_cid[0] if len(unique_cid) == 1 else unique_cid
                self.CAS = unique_cas if unique_cas else None
                self.InChi = unique_inchi[0] if unique_inchi else None
                self.InChiKey = unique_inchikey[0] if unique_inchikey else None
                self.M_array = arr_m

                # Select the record with the minimum M (if available)
                if np.isnan(arr_m).all():
                    self.M = None
                    self.formula = None
                    self.smiles = None
                else:
                    idx_min = np.nanargmin(arr_m)
                    self.M = arr_m[idx_min]
                    self.formula = all_data["formula"][idx_min]
                    self.smiles = all_data["smiles"][idx_min]

                # Process logP: store only valid (non-NaN) values if available
                valid_logp = arr_logp[~np.isnan(arr_logp)]
                self.logP = valid_logp if valid_logp.size > 0 else None

                # Store estimated on the van-der-Waals volume
                # (they are all wrong, then two estimates)
                self.vdWvolume = self.volume_3d # for future use (Dwelle use it)
                alpha = 1/0.65 if self.count_rings["aromatic"]>0 else 1.3
                self.vdWvolume2 = 1/alpha * self.molarvolumeMiller * 10.0/6.02214076 # in A3


        # Case (b): name is None, M is provided => generic substance
        # ----------------------------------------------------------------
        elif (name is None) and (M_array is not None):
            # No search performed
            if M_array.size == 1:
                self.compound = "single molecular weight"
            else:
                self.compound = (f"list of molecular weights ranging from "
                                 f"{float(np.min(M_array))} to {float(np.max(M_array))}")

            # name => "generic" or if user explicitly set name=..., handle it here
            self.name = "generic"  # from instructions
            self.cid = None
            self.CAS = None
            self.InChi = None
            self.InChiKey = None
            self.M_array = M_array
            self.M = float(np.min(M_array))
            self.formula = None
            self.smiles = None
            self.logP = logP_array  # user-supplied or None
            self.vdWvolume, self.vdWvolume2 = None, None

        # Case (c): name is not None and M is provided => surrogate
        # ----------------------------------------------------------------
        elif (name is not None) and (M_array is not None):
            # No search is done, it doesn't exist in PubChem
            if M_array.size == 1:
                self.compound = "single molecular weight"
            else:
                self.compound = (f"list of molecular weights ranging from "
                                 f"{float(np.min(M_array))} to {float(np.max(M_array))}")

            self.name = name
            self.cid
            self.CAS = None
            self.InChi = None
            self.InChiKey = None
            self.M_array = M_array
            self.M = float(np.min(M_array))
            self.formula = None
            self.smiles = None
            self.logP = logP_array
            self.vdWvolume, self.vdWvolume2 = None, None

        else:
            # If none of these scenarios apply, user gave incomplete or conflicting args
            raise ValueError("Invalid arguments. Provide either name for search (case a), "
                             "or M for a generic (case b), or both for a surrogate (case c).")


        # Model validation and paramameterization
        # ------------------------------------------------
        # migrant parameters are populated to the template
        # ------------------------------------------------

        # Diffusivity model
        if Dmodel is not None:
            self._validate_and_set_model("D",Dmodel,Dtemplate,
                                         {"M": self.M, "logP": self.logP, "Vvdw": self.volumeDwelle},
                                         MigrationPropertyModels,MigrationPropertyModel_validator)
        else:
            self.D = None
            self.Dtemplate = None

        # Henry-like model
        if kmodel is not None:
            self._validate_and_set_model("k",kmodel,ktemplate,
                                         {"Pi": self.polarityindex, "Vi": self.molarvolumeMiller},
                                         MigrationPropertyModels,MigrationPropertyModel_validator)
        else:
            self.k = None
            self.ktemplate = None


    def _validate_and_set_model(self, prop, model, template, update_params,PropertyModel,PropertyModelValidator):
        """
        Generic method for validating and setting a migration property model.

        Parameters:
            prop (str): The property identifier (e.g. "D" or "k").
            model (str or None): The model name.
            template (dict or None): A dictionary template to be updated.
            update_params (dict): Extra parameters to update the template.
            PropertyModel: MigrationPropertyModels from patankar.property
            PropertyModelValidator: MigrationPropertyModel_validator from patankar.property
        """
        if model is not None:
            if not isinstance(model, str):
                raise TypeError(f"{prop}model should be str not a {type(model).__name__}")
            if model not in PropertyModel[prop]:
                raise ValueError(f'The {prop} model "{model}" does not exist')
            model_class = PropertyModel[prop][model]
            if not PropertyModelValidator(model_class, model, prop):
                raise TypeError(f'The {prop} model "{model}" is corrupted')
            if template is None:
                template = {}
            if not isinstance(template, dict):
                raise TypeError(f"{prop}template should be a dict not a {type(template).__name__}")
            setattr(self, prop, model_class)
            temp_copy = template.copy()
            temp_copy.update(update_params)
            setattr(self, f"{prop}template", temp_copy)
        else:
            setattr(self, prop, None)
            setattr(self, f"{prop}template", None)

    # helper property to combine D and Dtemplate
    @property
    def Deval(self):
        """Return a callable function that evaluates D with updated parameters."""
        if self.D is None:
            return lambda **kwargs: None  # Return a function that always returns None
        def func(**kwargs):
            updated_template = dict(self.Dtemplate, **kwargs)
            return self.D.evaluate(**updated_template)
        return func

    # helper property to combine k and ktemplate
    @property
    def keval(self):
        """Return a callable function that evaluates k with updated parameters."""
        if self.k is None:
            return lambda **kwargs: None  # Return a function that always returns None
        def func(**kwargs):
            updated_template = dict(self.ktemplate, **kwargs)
            return self.k.evaluate(**updated_template)
        return func



    def __repr__(self):
        """Formatted string representation summarizing key attributes."""
        # Define header
        info = [f"<{self.__class__.__name__} object>"]
        # Collect attributes
        attributes = {
            "Compound": self.compound,
            "Name": self.name,
            "cid": self.cid,
            "CAS": self.CAS,
            "M (min)": self.M,
            "M_array": self.M_array if self.M_array is not None else "N/A",
            "formula": self.formula,
            "smiles": self.smiles if hasattr(self,"smiles") else "N/A",
            "InChiKey": self.InChiKey if hasattr(self,"InChiKey") else "N/A",
            "logP": self.logP,
            "P' (calc)": self.polarityindex
        }
        if isinstance(self,migrantToxtree):
            attributes["Compound"] = self.ToxTree["IUPACTraditionalName"]
            attributes["Name"] = self.ToxTree["IUPACName"]
            attributes["Toxicology"] = self.CramerClass
            attributes["TTC"] = f"{self.TTC} {self.TTCunits}"
            attributes["CF TTC"] = f"{self.CFTTC} {self.CFTTCunits}"
            alerts = self.alerts
            # Process alerts
            alert_index = 0
            for key, value in alerts.items():
                if key.startswith("Alert") and key != "Alertscounter" and value.upper() == "YES":
                    alert_index += 1
                    # Convert key name to readable format (split at capital letters)
                    alert_text = ''.join([' ' + char if char.isupper() and i > 0 else char for i, char in enumerate(key)])
                    attributes[f"alert {alert_index}"] = alert_text.strip()  # Remove leading space

        # Determine column width based on longest attribute name
        key_width = max(len(k) for k in attributes.keys()) + 2  # Add padding
        # Format attributes with indentation
        for key, value in attributes.items():
            formatted_key = f"{key}:".rjust(key_width)
            formatted_value = self.dispmax(value)
            info.append(f"  {formatted_key} {formatted_value}")
        # Print formatted representation
        repr_str = "\n".join(info)
        print(repr_str)
        # Return a short summary for interactive use
        return str(self)

    def __str__(self):
        """Formatted string representing the migrant"""
        onename = self.name[0] if isinstance(self.name,list) else self.name
        return f"<{self.__class__.__name__}: {self.dispmax(onename,16)} - M={self.M} g/mol>"

    def dispmax(self,content,maxwidth=None):
        """ optimize display """
        strcontent = str(content)
        maxwidth = self._maxdisplay if maxwidth is None else min(maxwidth,self._maxdisplay)
        if len(strcontent)>maxwidth:
            nchar = round(maxwidth/2)
            return strcontent[:nchar]+" [...] "+strcontent[-nchar:]
        else:
            return content

    # calculated propeties (rough estimates)
    @property
    def polarityindex(self,logP=None,V=None):
        """
            Computes the polarity index (P') of the compound.

            The polarity index (P') is derived from the compound's logP value and
            its molar volume V(), using an empirical (fitted) quadratic equation:

                E = logP * ln(10) - S
                P' = (-B - sqrt(B² - 4A(C - E))) / (2A)

            where:
                - S is the entropy contribution, calculated from molar volume.
                - A, B, C are empirical coefficients.

            Returns
            -------
            float
                The estimated polarity index P' based on logP and molar volume.

            Notes
            -----
            - For highly polar solvents (beyond water), P' saturates at **10.2**.
            - For extremely hydrophobic solvents (beyond n-Hexane), P' is **0**.
            - Accuracy is dependent on the reliability of logP and molar volume models.

            Example
            -------
            >>> compound.polarityindex
            8.34  # Example output
        """
        return polarity_index(logP=self.logP if logP is None else logP,
                              V=self.molarvolumeMiller if V is None else V)

    @property
    def molarvolumeMiller(self, a=0.997, b=1.03):
        """
        Estimates molar volume using the Miller empirical model.

        The molar volume (V_m) is calculated based on molecular weight (M)
        using the empirical formula:

            V_m = a * M^b  (cm³/mol)

        where:
            - `a = 0.997`, `b = 1.03` are empirically derived constants.
            - `M` is the molecular weight (g/mol).
            - `V_m` is the molar volume (cm³/mol).

        Returns
        -------
        float
            Estimated molar volume in cm³/mol.

        Notes
        -----
        - This is an approximate model and may not be accurate for all compounds.
        - Alternative models include the **Yalkowsky & Valvani method**.

        Example
        -------
        >>> compound.molarvolumeMiller
        130.5  # Example output
        """
        return a * self.M**b


    @property
    def molarvolumeLinear(self):
        """
        Estimates molar volume using a simple linear approximation.

        This method provides a rough estimate of molar volume, particularly
        useful for small to mid-sized non-ionic organic molecules. It is based on:

            V_m = 0.935 * M + 14.2  (cm³/mol)

        where:
            - `M` is the molecular weight (g/mol).
            - `V_m` is the estimated molar volume (cm³/mol).
            - Empirical coefficients are derived from **Yalkowsky & Valvani (1980s)**.

        Returns
        -------
        float
            Estimated molar volume in cm³/mol.

        Notes
        -----
        - This method is often *"okay"* for non-ionic organic compounds.
        - Accuracy decreases for very large, ionic, or highly branched molecules.
        - More precise alternatives include **Miller's model** or **group contribution methods**.

        Example
        -------
        >>> compound.molarvolumeLinear
        120.7  # Example output
        """
        return 0.935 * self.M + 14.2

    # suggest an alternative D model
    def suggest_alt_Dmodel(self,material=None,index=0,RaiseError=True,RaiseWarning=True,**template):
        """suggest an alternative Dmodel based on Dmodel_extensions"""
        # local dependencies
        from patankar.property import PropertyModelSelector
        # special case when material is a str
        if isinstance(material,str):
            import patankar.layer as allmaterials
            if not hasattr(allmaterials,material):
                raise ValueError(f"{material} is not a valid class of patankar.layer")
            cls = getattr(allmaterials,material)
            material = cls()
        if RaiseError: # check all args
            from patankar.layer import layer
            if material is None:
                raise ValueError("material must be provided")
            if not isinstance(material, layer):
                # note that forcing layer.__compute_Dmodel(...RaiseError=True) may raise errors
                # since a different reference to layer might be in memory
                raise TypeError(f"material must be a layer not a {type(material).__name__}")
            if not(isinstance(index,int)):
                raise TypeError(f"index must be int not a {type(index).__name__}")
            if index>len(material):
                raise ValueError(f"index value {index} exceeds the number of layers {len(material)}")
        # build objects on which rules will be tested
        objects = {"material":material,"migrant":self,"index":index}
        # check Dmodels with objects
        availableExtensions = list(Dmodel_extensions.keys())
        applicableExtensions = [False]*len(availableExtensions)
        import patankar.property as allproperties
        for imodel,modelCode in enumerate(availableExtensions):
            modelObjects = [objects[o] for o in Dmodel_extensions[modelCode]["objects"]]
            modelRules = Dmodel_extensions[modelCode]["rules"].copy()
            modelRules_layer = modelRules[0]["list"][1] # 0=polymer rules, 1=polymer type
            modelRules_layer["index"] = objects["index"] # layer index
            applicableExtensions[imodel] = PropertyModelSelector(modelRules,modelObjects)

            # if the model is OK we try to import it and check at temperature T
            if applicableExtensions[imodel]:
                try:
                    modelclass = getattr(allproperties, availableExtensions[imodel])
                except AttributeError:
                    if RaiseWarning:
                        print(f"WARNING:: The D model {modelclass} is not defined in patankar.property.")
                    applicableExtensions[imodel] = False # the model is not avaiable
                    continue
                modeltest = modelclass.evaluate(**template) # we call the static method without instantiation
                if modeltest is None:
                    applicableExtensions[imodel] = False # the model is not suitable (we will take the next)

        first_true_index = next((i for i, val in enumerate(applicableExtensions) if val), None)
        # returns the name of the first applicable model, if not None
        return availableExtensions[first_true_index] if first_true_index is not None else None

    # suggest an alternative D class
    def suggest_alt_Dclass(self,material=None,index=0,RaiseError=True,RaiseWarning=True,**template):
        """returns an alternative Dclass based on Dmodel_extensions"""
        alt_classname = self.suggest_alt_Dmodel(material=material,index=index,RaiseError=RaiseError,RaiseWarning=RaiseWarning,**template)
        if alt_classname is None:
            return None
        import patankar.property as allproperties
        try:
            alt_class = getattr(allproperties, alt_classname)
        except AttributeError:
            if RaiseWarning: # it is a double check
                print(f"WARNING:: The D model {alt_classname} is not defined in patankar.property.")
            return None
        if not isinstance(alt_class, type):
            raise TypeError(f"Expected a class for {alt_classname}, but found {type(alt_class).__name__}.")
        return alt_class

    # suggest an alternative D class
    def check_alt_propclass(self,alt_classname):
        """returns True if a class property exists in patankar.property"""
        if alt_classname is None:
            return False
        import patankar.property as allproperties
        try:
            alt_class = getattr(allproperties, alt_classname)
        except AttributeError:
            return False
        if not isinstance(alt_class, type):
            return False
        return True

    # --------------------------------------------------------------------------
    #            [   M O L E C U L A R   D E S C R I P T O R S   ]
    #         experimental implementation to remove external dependencies
    #                               results are AS IS
    # --------------------------------------------------------------------------
    @property
    def count_rings(self):
        """
        Count aromatic and non-aromatic rings separately from canonical SMILES.

        The method removes bracketed expressions to avoid counting digits
        that are part of atomic specifications, then finds ring closure digits.
        It handles both single-digit ring closures and multi-digit closures (e.g., %12).

        The method is very rough, do not expect the best results.

        Returns:
            dict: {
                'total': int,
                'aromatic': int,
                'non_aromatic': int,
                'fusions': {
                    'AlAr': int,
                    'AlAl': int,
                    'ArAr': int,
                    'S': int
                }
            }
        """
        smiles = self.smiles
        smiles_no_brackets = re.sub(r'\[.*?\]', '', smiles)

        closures = {}
        for match in re.finditer(r'(\%\d{2}|\d)', smiles_no_brackets):
            marker = match.group()
            pos = match.start()
            closures.setdefault(marker, []).append(pos)

        aromatic_rings = set()
        non_aromatic_rings = set()
        ring_aromaticity = {}

        for marker, positions in closures.items():
            if len(positions) == 2:
                start, end = sorted(positions)
                ring_substring = smiles_no_brackets[start:end+len(marker)]

                is_aromatic = (
                    bool(re.search(r'[a-z]', ring_substring)) or
                    len(re.findall('=', ring_substring)) >= 2
                )

                ring_aromaticity[marker] = 'Ar' if is_aromatic else 'Al'

                if is_aromatic:
                    aromatic_rings.add(marker)
                else:
                    non_aromatic_rings.add(marker)

        fusion_counts = {'AlAl': 0, 'AlAr': 0, 'ArAr': 0, 'S': 0}
        markers = list(closures.keys())

        for i in range(len(markers)):
            for j in range(i + 1, len(markers)):
                # Rings are fused if they share a closure digit (atom)
                if len(set(closures[markers[i]]) & set(closures[markers[j]])) >= 1:
                    types = sorted([ring_aromaticity[markers[i]], ring_aromaticity[markers[j]]])
                    fusion_type = ''.join(types)
                    if fusion_type in fusion_counts:
                        fusion_counts[fusion_type] += 1

        # Special pattern (S-type) detection simplified
        fusion_counts['S'] = len(re.findall(r's1.*?c.*?c.*?c.*?c.*?c.*?c1.*?c', smiles_no_brackets))

        return {
            'total': len(aromatic_rings) + len(non_aromatic_rings),
            'aromatic': len(aromatic_rings),
            'non_aromatic': len(non_aromatic_rings),
            'fusions': fusion_counts
        }



    @property
    def parse_formula(self):
        """
        Parse the molecular formula into a dictionary of element counts.

        For example, 'C15H16O2' will be parsed as:
          {'C': 15, 'H': 16, 'O': 2}

        Returns:
          dict: Dictionary mapping element symbols to their counts.
        """
        formula = self.formula
        pattern = re.compile(r'([A-Z][a-z]*)(\d*)')
        counts = {}
        for (element, count) in pattern.findall(formula):
            count = int(count) if count else 1
            counts[element] = counts.get(element, 0) + count
        return counts

    @property
    def volume_3d(self):
        """
        Compute the molecular 3D van-der-Waals volume using a linear additive model.

        This function applies an empirical, linear-additive scheme:
        Adapted code from: https://github.com/hachmannlab/Slonimskii_method_vdW/blob/master/Slonimskii_method_vdW.py. The intent is to remove the dependence to babel.

        \[
        \begin{aligned}
        V_{vdw} =\,& n_H \times 7.24 + n_C \times 20.58 + n_N \times 15.6 + n_O \times 14.71 + n_F \times 13.31 \\
        &+ n_{Cl} \times 22.45 + n_{Br} \times 26.52 + n_{I} \times 32.52 + n_P \times 24.43 \\
        &+ n_S \times 24.43 + n_{Si} \times 38.79 + n_{Se} \times 28.73 + n_{Te} \times 36.62 \\
        &- 5.92 \times (\text{bonds}) - 14.7 \times (\text{no\_ar}) - 3.8 \times (\text{no\_non\_ar}) \\
        &+ 5 \times (\text{no\_f\_ring\_ArAr}) + 3 \times (\text{no\_f\_ring\_AlAr}) + 1 \times (\text{no\_f\_ring\_AlAl}) \\
        &- 5 \times (\text{no\_f\_ring\_S})
        \end{aligned}
        \]

        Here, the bond count is estimated using the heuristic:

        \[
        \text{bonds} \approx (\text{number of heavy atoms} - 1) + (\text{number of rings})
        \]

        and heavy atoms are all atoms except hydrogen.

        For this implementation:
          - All rings detected in the SMILES are assumed to be aromatic (valid for bisphenol A).
          - Fused ring corrections (no_f_ring_*) are set to zero.

        Parameters:
          smiles (str): Canonical SMILES string of the molecule.
          formula (str): Molecular formula (e.g., 'C15H16O2') to retrieve accurate hydrogen counts.

        Returns:
          float: Estimated 3D van-der-Waals volume in Å³.
        """
        # Parse the molecular formula to obtain element counts.
        counts = self.parse_formula
        no_H = counts.get('H', 0)
        no_C = counts.get('C', 0)
        no_N = counts.get('N', 0)
        no_O = counts.get('O', 0)
        no_F = counts.get('F', 0)
        no_Cl = counts.get('Cl', 0)
        no_Br = counts.get('Br', 0)
        no_I = counts.get('I', 0)
        no_P = counts.get('P', 0)
        no_S = counts.get('S', 0)
        no_Si = counts.get('Si', 0)
        no_Se = counts.get('Se', 0)
        no_Te = counts.get('Te', 0)

        # Count rings using the SMILES string.
        num_rings = self.count_rings

        # Assume all rings are aromatic (for bisphenol A, both rings are aromatic).
        no_ar = num_rings["aromatic"]
        no_non_ar = num_rings["non_aromatic"]
        no_total = num_rings["total"]

        # Estimate the number of bonds.
        # Heavy atoms are those except hydrogen.
        heavy_atoms = no_C + no_N + no_O + no_F + no_Cl + no_Br + no_I + no_P + no_S + no_Si + no_Se + no_Te
        # For a tree structure, bonds = (heavy_atoms - 1). Each ring adds one extra bond.
        bonds = (heavy_atoms - 1) + no_total

        # For fused ring corrections
        no_f_ring_ArAr = num_rings["fusions"]["ArAr"]
        no_f_ring_AlAr = num_rings["fusions"]["AlAr"]
        no_f_ring_AlAl = num_rings["fusions"]["AlAl"]
        no_f_ring_S = num_rings["fusions"]["S"]
        # remove 1 for each fused ring
        bonds = bonds - no_f_ring_ArAr - no_f_ring_AlAr - no_f_ring_AlAl - no_f_ring_S

        # Compute the initial volume from atomic contributions.
        V_vdw = (no_H)*7.24 + (no_C)*20.58 + (no_N)*15.6 + (no_O)*14.71 + (no_F)*13.31 + \
                (no_Cl)*22.45 + (no_Br)*26.52 + (no_I)*32.52 + (no_P)*24.43 + (no_S)*24.43 + \
                (no_Si)*38.79 + (no_Se)*28.73 + (no_Te)*36.62

        # Apply corrections for bonds and rings.
        V_vdw = V_vdw - 5.92*(bonds) - 14.7*(no_ar) - 3.8*(no_non_ar) + \
                5*(no_f_ring_ArAr) + 3*(no_f_ring_AlAr) + 1*(no_f_ring_AlAl) - 5*(no_f_ring_S)

        return V_vdw

    @property
    def volumeDwelle(self):
        """
            Returns the approximate volume as estimated by Dwelle
            In the original paper [1], molecular volumes are calculated with [2].
            The values reported in [1] are particularly low such as the H contribution
            was missing. As a result a linear correction is proposed based on phenanthrene.

            [1] https://onlinelibrary.wiley.com/doi/full/10.1002/pts.2638
            [2] https://www.molinspiration.com/services/volume.html

            The minimum of vdWvolume (calculated with volume_3D) and vdWvolume2 (inferred
            from molarvolumeMiller using a linear correlation for aromatic and non-aromatic
            molecules) is taken. This procedure aims at preserving conservatism in Dwelle
            estimates.

        """
        return 172.2/221.1 * min(self.vdWvolume,self.vdWvolume2) # 172.2/221.7 = 0.7767

# %% Class migrantToxtree extending migrant class with Toxtree data
"""
===============================================================================
SFPPy loadpubchem extension: Interface to Toxtree
===============================================================================
This module extends SFPPy migrants with the capability to interface with Toxtree for toxicological assessments:
- Cramer classification
- Detection of structural alerts from SMILES

The module leverages existing PubChem data managed by SFPPy.

===============================================================================
"""
class migrantToxtree(migrant):
    """
    Extends the `migrant` class to integrate Toxtree for toxicological assessments.
    This class retrieves chemical data from PubChem, caches results, and runs Toxtree
    for toxicological classification.

    Features:
    - Downloads and caches molecular structure files (SDF) from PubChem.
    - Interfaces with Toxtree to perform Cramer classification and detect toxicological alerts.
    - Implements a multi-level cache system for efficiency.
    - Cleans and standardizes field names for output consistency.
    - Provides control flags for cache refresh and regeneration.

    Attributes:
        PUBCHEM_ROOT_URL (str): Base URL for PubChem compound data.
        IMAGE_SIZE (tuple): Default image size for structure images.
        TOXTREE_ENGINES (dict): Mapping of Toxtree engines to class names.
        TOX_CLASSIFICATION (dict): Mapping of classification engines.
        cache_folder (str): Directory to store cached Toxtree results.
        structure_folder (str): Directory to store cached structure files.
        structure_file (str): Path to the SDF file for the compound.
        image_file (str): Path to the PNG image of the compound.
        refresh (bool): If True, forces Toxtree reprocessing from CSV.
        no_cache (bool): If True, forces full cache regeneration.

    Methods:
        __init__(compound_name, cache_folder='cache.ToxTree', structure_folder='structure', refresh=False, no_cache=False)
            Initializes the class, retrieves chemical data, and runs Toxtree default classification.

        _download_pubchem_data()
            Downloads and caches the SDF structure file from PubChem.

        _clean_field_names(data)
            Cleans field names by removing 'PUBCHEM_' prefix and applying CamelCase.

        _run_toxtree(engine)
            Runs Toxtree for the specified engine, handling caching and errors.

    Properties:
        cramer: Runs Toxtree with the Cramer classification engine.
        cramer2: Runs Toxtree with the Cramer2 classification engine.
        cramer3: Runs Toxtree with the Cramer3 classification engine.
        alerts: Runs Toxtree to detect toxicological alerts.
        has_alerts: Checks if any toxicological alerts were detected.

    Example:
        >>> substance = migrantToxtree("limonene")
        >>> c = substance.cramer
        >>> c2 = substance.cramer2
        >>> c3 = substance.cramer3
        >>> print("Cramer Class:", c)
        >>> print("Cramer2 Class:", c2)
        >>> print("Cramer3 Class:", c3)
    """

    PUBCHEM_ROOT_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
    IMAGE_SIZE = (640, 480)

    TOXTREE_ENGINES = {
        "default": "",
        "cramer": "toxTree.tree.cramer.CramerRules",
        "cramer2": "cramer2.CramerRulesWithExtensions",
        "cramer3": "toxtree.tree.cramer3.RevisedCramerDecisionTree",
        "kroes": "toxtree.plugins.kroes.Kroes1Tree",
        "dnabinding": "toxtree.plugins.dnabinding.DNABindingPlugin",
        "skin": "toxtree.plugins.skinsensitisation.SkinSensitisationPlugin",
        "eye": "eye.EyeIrritationRules",
        "Ames": "toxtree.plugins.ames.AmesMutagenicityRules"
    }

    TOX_CLASSIFICATION = {
        "default": "CramerRules",
        "cramer": "CramerRules",
        "cramer2": "CramerRules_WithExtensions",
        "cramer3": "RevisedCDT",
        "kroes": "KroesTTCDecisionTree"
    }

    TTC = [0.0025, 1.5, 9.0 , 30] # µg/kg bw/day
    TTCunits = "[µg/kg bw/day]"
    CFTTC = [ttc * 60 * 1 * 1e-3 for ttc in TTC] # mg/kg intake
    CFTTCunits = "[mg/kg food intake]"

    def __init__(self, compound_name, cache_folder='cache.ToxTree', structure_folder='structure', refresh=False, no_cache=False):
        super().__init__(compound_name)

        if isinstance(self.cid, list):
            if len(self.cid) != 1:
                raise ValueError(f"Multiple CIDs found for {compound_name}. Provide a unique compound.")
            self.cid = self.cid[0]
        if isinstance(self.smiles, list):
            if len(self.smiles) != 1:
                raise ValueError(f"Multiple SMILES found for {compound_name}. Provide a unique SMILES.")
            self.smiles = self.smiles[0]

        self.toxtree_root = os.path.join(_PATANKAR_FOLDER, 'private', 'toxtree')
        self.jar_path = os.path.join(self.toxtree_root, 'Toxtree-3.1.0.1851.jar')

        if not os.path.isfile(self.jar_path):
            raise FileNotFoundError(
                f"The Toxtree executable '{self.jar_path}' cannot be found.\n"
                f"Please follow the instructions in the README.md file located at '{self.toxtree_root}'."
            )

        self.cache_folder = os.path.join(_PATANKAR_FOLDER, cache_folder)
        self.structure_folder = os.path.join(self.cache_folder, structure_folder)
        os.makedirs(self.cache_folder, exist_ok=True)
        os.makedirs(self.structure_folder, exist_ok=True)

        self.structure_file = os.path.join(self.structure_folder, f'{self.cid}.sdf')
        self.image_file = os.path.join(self.structure_folder, f'{self.cid}.png')
        self.refresh = refresh
        self.no_cache = no_cache

        self._download_pubchem_data()
        tmp = self._run_toxtree("default")
        tmp["CramerValue"]=self.class_roman_to_int(tmp["CramerRules"])
        self.ToxTree = tmp
        self.CramerValue = tmp["CramerValue"]
        self.CramerClass = tmp["CramerRules"]
        self.TTC = self.TTC[self.CramerValue]
        self.CFTTC = self.CFTTC[self.CramerValue]

    def _download_pubchem_data(self):
        """Downloads and caches the SDF structure file and PNG thumbnail from PubChem."""
        if not os.path.isfile(self.structure_file) or self.no_cache:
            sdf_url = f"{self.PUBCHEM_ROOT_URL}/CID/{self.cid}/SDF"
            response = requests.get(sdf_url, timeout=60)
            if response.status_code == 200:
                with open(self.structure_file, 'wb') as f:
                    f.write(response.content)
            else:
                raise ValueError(f"Failed to download SDF file for CID {self.cid}.")
        if not os.path.isfile(self.image_file) or self.no_cache:
            png_url = f"{self.PUBCHEM_ROOT_URL}/CID/{self.cid}/PNG?image_size={self.IMAGE_SIZE[0]}x{self.IMAGE_SIZE[1]}"
            response = requests.get(png_url, timeout=60)
            if response.status_code == 200:
                with open(self.image_file, 'wb') as f:
                    f.write(response.content)
                self._crop_image(self.image_file)

    def _crop_image(self, image_path):
        """Crops white background from the PNG image."""
        if not PIL_AVAILABLE:
            return
        img = Image.open(image_path)
        img = img.convert("RGBA")
        bbox = img.getbbox()
        if bbox:
            img = img.crop(bbox)
            img.save(image_path)


    def _clean_field_names(self, data):
        """Cleans field names by removing PUBCHEM_, splitting with multiple delimiters, and capitalizing each word."""
        cleaned_data = {}
        substitutions = {
            "iupac": "IUPAC",
            "logp": "logP",
            "cramer": "Cramer",
            "cid": "CID",
            "inchi": "InChi",
            }
        # Create a case-insensitive regex pattern
        pattern = re.compile("|".join(re.escape(k) for k in substitutions.keys()), re.IGNORECASE)
        # Function for case-insensitive substitution
        def replace_case_insensitive(match):
            return substitutions[match.group(0).lower()]  # Lookup in lowercase, replace as defined
        for key, value in data.items():
            # Remove PUBCHEM_ prefix
            if key.startswith("PUBCHEM_"):
                key = key.replace("PUBCHEM_", "")
            # Split using multiple delimiters (space, underscore, comma)
            words = re.split(r'[ ,_:\.]+', key)
            # Capitalize each word
            cleaned_key = ''.join(word.capitalize() for word in words)
            # Apply case-sensitive substitutions
            cleaned_key = pattern.sub(replace_case_insensitive, cleaned_key)
            # Store in the cleaned dictionary
            cleaned_data[cleaned_key] = value
        return cleaned_data

    def _run_toxtree(self, engine):
        if engine not in self.TOXTREE_ENGINES:
            raise ValueError(f"Unknown Toxtree engine: {engine}")

        engine_class = self.TOXTREE_ENGINES[engine]
        csv_file = os.path.join(self.cache_folder, f'{self.cid}.{engine}.csv')
        json_file = os.path.join(self.cache_folder, f'{self.cid}.{engine}.json')

        if os.path.isfile(json_file) and not self.refresh:
            with open(json_file, 'r') as f:
                return json.load(f)

        if not os.path.isfile(csv_file) or self.no_cache:
            cmd = ['java', '-jar', self.jar_path, '-n', '-i', self.structure_file, '-o', csv_file]
            if engine_class:
                cmd.extend(['-m', engine_class])
            try:
                current_dir = os.getcwd()
                os.chdir(self.toxtree_root) # toxtree needs to run from its installation folder
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                os.chdir(current_dir) # restore path
                if not os.path.isfile(csv_file):
                    raise RuntimeError(
                        f"Error: Toxtree failed to generate the output file for {engine}.\n"
                        f"Command: {' '.join(cmd)}\n"
                        f"Output: {result.stdout}\n"
                        f"Error: {result.stderr}"
                    )
            except subprocess.CalledProcessError as e:
                os.chdir(current_dir) # restore path
                raise RuntimeError(
                    f"Error executing Toxtree for {engine}.\n"
                    f"Command: {' '.join(cmd)}\n"
                    f"Output: {e.stdout}\n"
                    f"Error: {e.stderr}"
                )

        df = pd.read_csv(csv_file)
        cleaned_data = self._clean_field_names(df.to_dict(orient='records')[0]) if not df.empty else {}
        with open(json_file, 'w') as f:
            json.dump(cleaned_data, f, indent=4)
        return cleaned_data


    def class_roman_to_int(self,text):
        """Converts 'Class X' (where X is I, II, III, IV, V) into an integer (1-5), case insensitive."""
        # Define mapping of Roman numerals to integers
        roman_to_int = {
            "I": 1, "II": 2, "III": 3, "IV": 4, "V": 5
        }
        # Regex pattern to detect 'Class X' with valid Roman numerals
        pattern = re.compile(r'\bClass\s+(I{1,3}|IV|V)\b', re.IGNORECASE)
        # Search for pattern in text
        match = pattern.search(text)
        if match:
            roman_numeral = match.group(1).upper()  # Extract and normalize the numeral
            return roman_to_int.get(roman_numeral, None)  # Convert to int
        return None  # Return None if no valid match


    @property
    def cramer(self):
        tmp = self._run_toxtree('cramer')
        tmp["CramerValue"]=self.class_roman_to_int(tmp["CramerRules"])
        tmp["IsCramerFlag"] = not pd.isna(tmp["Cramerflags"]) if "Cramerflags" in tmp else None
        return tmp

    @property
    def cramer2(self):
        tmp = self._run_toxtree('cramer2')
        tmp["CramerValue"] = self.class_roman_to_int(tmp["CramerRulesWithExtensions"]) \
            if "CramerRulesWithExtensions" in tmp else None
        tmp["IsCramerFlag"] = not pd.isna(tmp["Cramerflags"]) if "Cramerflags" in tmp else None
        return tmp

    @property
    def cramer3(self):
        tmp = self._run_toxtree('cramer3')
        tmp["CramerValue"] = self.class_roman_to_int(tmp["CramerRulesWithExtensions"]) \
            if "CramerRulesWithExtensions" in tmp else None
        tmp["IsCramerFlag"] = not pd.isna(tmp["Cramerflags"]) if "Cramerflags" in tmp else None
        return tmp

    @property
    def alerts(self):
        return self._run_toxtree('skin')

    @property
    def has_alerts(self):
        return len(self.alerts) > 0



# %% debug
# ==========================
# Usage example:
# ==========================
if __name__ == "__main__":
    # debug
    m=migrant("bisphenol A")
    m.count_rings
    m.volume_3d
    m=migrant("water")
    m.polarityindex
    # examples
    db = CompoundIndex()
    df_simple = db.find("limonene", output_format="simple")
    df_simple = db.find("aspirin", output_format="simple")
    df_simple = db.find("irganox 1076", output_format="simple")
    df_simple = db.find("anisole", output_format="simple")
    print("Simple result:\n", df_simple)

    df_full = db.find("anisole", output_format="full")
    print("Full result:\n", df_full)

    # for migration modeling
    m = migrant(name='anisole')
    print(m)
    m = migrant(name='limonene')
    print(m)
    m = migrant(name='irganox 1076')
    print(m)
    m = migrant(name='irgafos 168')
    print(m)
    m = migrant() # toluene
    print(m)
    # Piringer D value (several models can be implemented in module property.py)
    Dval = m.Deval(polymer="PET",T=20)
    print(Dval)

    # MigranToxtree tests
    substance = migrantToxtree("irganox 1010")
    c = substance.cramer
    c2 = substance.cramer2
    c3 = substance.cramer3
    print("Cramer Class:", c)
    print("Cramer2 Class:", c2)
    print("Cramer3 Class:", c3)

    # suggest an alternative D model
    from patankar.layer import gPET, LDPE, PP, rigidPVC
    material = gPET()+LDPE()+PP()+rigidPVC()
    m.suggest_alt_Dmodel(material,3)
    m.suggest_alt_Dmodel(material,1)
