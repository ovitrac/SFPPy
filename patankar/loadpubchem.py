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


@version: 1.2
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2024-02-17
@rev: 2025-02-21

Version History
---------------
- 1.0: Initial version, supporting local caching, synonyms index, and direct PubChem lookup.

"""


import os
import json
import re
import glob
import pandas as pd
import numpy as np
from datetime import datetime

# private version of pubchempy
from private.pubchempy import get_compounds

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
        self.cache_dir = cache_dir
        os.makedirs(self.cache_dir, exist_ok=True)

        self.index_file = os.path.join(cache_dir, index_file)
        # Regex to identify CAS-like strings, e.g. "1234-56-7"
        self._cas_regex = re.compile(r'^\d{1,7}-\d{2}-\d$')

        # Attempt to load existing index; if missing or invalid, rebuild
        if not os.path.isfile(self.index_file):
            self.refresh_index()
        else:
            with open(self.index_file, "r", encoding="utf-8") as f:
                try:
                    self.index = json.load(f)
                except json.JSONDecodeError:
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
        qlower = query.strip().lower()

        if qlower not in self.index:
            # Not found locally => do a PubChem call
            matches = get_compounds(query, 'name')
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

# %% Class migrant (high-level)

# Main compound database
dbdefault = CompoundIndex(cache_dir="cache.PubChem", index_file="pubchem_index.json")

# Migrant class
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

    def __init__(self, name=None, M=None, logP=None,
                 Dmodel = "Piringer",
                 Dtemplate = {"polymer":"LLDPE", "M":50, "T":40}, # do not use None
                 kmodel = None,
                 ktemplate = {},
                 db=dbdefault):
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
            - Set a Henri-like model using
                    - kmodel="model name"
                      default =None
                    - ktemplate=template dict coding for the key:value parameters
                      default = {}
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
        from property import MigrationPropertyModels, MigrationPropertyModel_validator

        self.compound = None   # str
        self.name = None       # str or list
        self.CAS = None        # list or None
        self.M = None          # float
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
                # No record found
                self.compound = name
                self.name = [name]
                self.CAS = []
                self.M_array = np.array([], dtype=float)
                self.M = None
                self.formula = None
                self.logP = None
            else:
                # Possibly multiple matching rows
                self.compound = name

                all_names = []
                all_cas = []
                all_m = []
                all_formulas = []
                all_logP = []

                for _, row in df.iterrows():
                    # Gather a list/set of names
                    row_names = row.get("name", [])
                    if isinstance(row_names, str):
                        row_names = [row_names]
                    row_syns = row.get("synonyms", [])
                    combined_names = set(row_names) | set(row_syns)
                    all_names.extend(list(combined_names))

                    # CAS
                    row_cas = row.get("CAS", [])
                    if row_cas:
                        all_cas.extend(row_cas)

                    # M
                    row_m = row.get("M", None)
                    if row_m is not None:
                        try:
                            all_m.append(float(row_m))
                        except:
                            all_m.append(np.nan)
                    else:
                        all_m.append(np.nan)

                    # logP
                    row_logp = row.get("logP", None)
                    if row_logp not in (None, ""):
                        try:
                            all_logP.append(float(row_logp))
                        except:
                            all_logP.append(np.nan)
                    else:
                        all_logP.append(np.nan)

                    # formula (as a string)
                    row_formula = row.get("formula", None)
                    # Even if None, we append so the index lines up with M
                    all_formulas.append(row_formula)

                # Convert to arrays
                arr_m = np.array(all_m, dtype=float)
                arr_logp = np.array(all_logP, dtype=float)

                # Some dedup / cleaning
                unique_names = list(set(all_names))
                unique_cas = list(set(all_cas))

                # Store results in the migrant object
                self.name = unique_names
                self.CAS = unique_cas if unique_cas else None
                self.M_array = arr_m
                # Minimum M
                if np.isnan(arr_m).all():
                    self.M = None
                    self.formula = None
                else:
                    idx_min = np.nanargmin(arr_m)    # index of min M
                    self.M = arr_m[idx_min]         # pick that M
                    self.formula = all_formulas[idx_min]  # pick formula from same record

                # Valid logP
                valid_logp = arr_logp[~np.isnan(arr_logp)]
                if valid_logp.size > 0:
                    self.logP = valid_logp  # or store as a list/mean/etc.
                else:
                    self.logP = None

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
            self.CAS = None
            self.M_array = M_array
            self.M = float(np.min(M_array))
            self.formula = None
            self.logP = logP_array  # user-supplied or None

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
            self.CAS = None
            self.M_array = M_array
            self.M = float(np.min(M_array))
            self.formula = None
            self.logP = logP_array

        else:
            # If none of these scenarios apply, user gave incomplete or conflicting args
            raise ValueError("Invalid arguments. Provide either name for search (case a), "
                             "or M for a generic (case b), or both for a surrogate (case c).")


        # Model validation and paramameterization
        # ----------------------------------------

        # Diffusivity model
        if Dmodel is not None:
            if not isinstance(Dmodel,str):
                raise TypeError(f"Dmodel should be str not a {type(Dmodel).__name__}")
            if Dmodel not in MigrationPropertyModels["D"]:
                raise ValueError(f'The diffusivity model "{Dmodel}" does not exist')
            Dmodelclass = MigrationPropertyModels["D"][Dmodel]
            if not MigrationPropertyModel_validator(Dmodelclass,Dmodel,"D"):
                raise TypeError(f'The diffusivity model "{Dmodel}" is corrupted')
            if Dtemplate is None:
                Dtemplate = {}
            if not isinstance(Dtemplate,dict):
                raise TypeError(f"Dtemplate should be a dict not a {type(Dtemplate).__name__}")
            self.D  = Dmodelclass
            self.Dtemplate = Dtemplate.copy()
            self.Dtemplate.update(M=self.M,logP=self.logP)
        else:
            self.D = None
            self.Dtemplate = None

        # Henri-like model
        if kmodel is not None:
            if not isinstance(kmodel,str):
                raise TypeError(f"kmodel should be str not a {type(kmodel).__name__}")
            if kmodel not in MigrationPropertyModels["k"]:
                raise ValueError(f'The Henri-like model "{kmodel}" does not exist')
            kmodelclass = MigrationPropertyModels["k"][kmodel]
            if not MigrationPropertyModel_validator(kmodelclass,Dmodel,"k"):
                raise TypeError(f'The Henri-like model "{kmodel}" is corrupted')
            if ktemplate is None:
                ktemplate = {}
            if not isinstance(ktemplate,dict):
                raise TypeError(f"ktemplate should be a dict not a {type(ktemplate).__name__}")
            self.k  = kmodelclass
            self.ktemplate = ktemplate.copy()
            self.ktemplate.update(M=self.M,logP=self.logP)
        else:
            self.k = None
            self.ktemplate = None

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
            "CAS": self.CAS,
            "M (min)": self.M,
            "M_array": self.M_array if self.M_array is not None else "N/A",
            "formula": self.formula,
            "logP": self.logP
        }
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


# ==========================
# Usage example:
# ==========================
if __name__ == "__main__":
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
