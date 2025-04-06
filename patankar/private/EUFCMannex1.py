"""
Module: eu_fcm_annex1

This module implements a robust database manager for Annex I of EU Regulation 10/2011 EC, which lists authorized food contact substances.
The module reads a CSV file (downloaded from ECHA) and converts each row into a JSON record stored in a cache directory (as files named recXXXXX.annex1.json).
It then builds a global index (stored as annex1_index.json in the cache folder) mapping key fields—such as chemical name, EC number, CAS number, FCM number, reference,
and PubChem CID—to their corresponding record numbers.

Key features include:
  - Conversion of CSV rows into JSON records with ordered keys.
  - Use of a primary key (the record number, which is 1-based and stored in the "record" field) instead of a zero-based Python index.
  - Lookup methods by CAS, record number, PubChem CID, as well as by other fields (byname, byEC, byFCM, byRef).
  - Additional search functions that filter records by migration limit (SML and SMLT) ranges.
  - Caching of substances with missing PubChem results in a separate file ("missing.pubchem.annex1.json") so that lookup requests for substances with a nonempty CAS that fail to resolve are not repeated.
  - Custom __repr__ and __str__ methods for the record objects, with text wrapping that allows the first line to wrap at a set width and subsequent lines to wrap at a separate width.
  - Enhanced error handling in __getitem__ that returns informative messages when the provided key (either record number or CAS) is out of range.

Global Rate Limiting:
  The module assumes that the PubChem query throttling is managed by a module-level global variable (e.g. PubChem_lastQueryTime) within the imported
  'migrant' function from the patankar.loadpubchem module. This state is maintained across multiple instances of EuFCMannex1.

Usage Example:
  >>> db = EuFCMannex1()
  >>> rec = db[1060]             # Retrieves the record with record number 1060.
  >>> rec_by_cas = db.byCAS("102-39-6")
  >>> for r in db:
  ...     print(r["name"])
  >>> if 10 in db:
  ...     print("Record 10 exists.")
  >>> if 113194 in db:
  ...     print("PubChem cid 113194 exists.")

@version: 1.40
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2024-01-10
@rev: 2025-03-27

"""

import os,csv, json, datetime, time, re, textwrap

__all__ = ['EuFCMannex1', 'annex1record', 'annex1record_ext', 'custom_wrap', 'printWARN']


__project__ = "SFPPy"
__author__ = "Olivier Vitrac"
__copyright__ = "Copyright 2022"
__credits__ = ["Olivier Vitrac"]
__license__ = "MIT"
__maintainer__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@agroparistech.fr"
__version__ = "1.41"

# Default value for SML when field is empty.
SMLdefault = 60.0

# Module-level variables to track last warning message and its timestamp
_LAST_WARN_ = None
_T_LAST_WARN_ = 0.0

# %% private functions
def custom_wrap(text, width=40, indent=" " * 22):
    # Wrap first line to first_width characters.
    first_line = textwrap.wrap(text, width=width)
    if not first_line:
        return ""
    # Get the first line and then compute the remaining text.
    first = first_line[0]
    remaining = text[len(first):].lstrip()
    # Wrap the remaining text to subsequent_width characters.
    subsequent_lines = textwrap.wrap(remaining, width=width)
    # Prepend the indent to subsequent lines.
    wrapped = [first] + [indent + line for line in subsequent_lines]
    return "\n".join(wrapped)

def printWARN(message: str, tsilent: float = 10.0):
    """
    Print a warning message only if:
    - it's different from the last one, or
    - more than `tsilent` seconds have passed since the last identical warning.

    Parameters:
    ----------
    message : str
        The warning message to display.
    tsilent : float, optional
        Minimum time (in seconds) between repeated identical warnings.
    """
    global _LAST_WARN_, _T_LAST_WARN_
    tnow = time.time()
    if message != _LAST_WARN_ or (tnow - _T_LAST_WARN_ > tsilent):
        print(message)
        _LAST_WARN_ = message
        _T_LAST_WARN_ = tnow

# %% main classes
# ----------------------------------------------------------------------
# annex1record: a subclass of dict to represent one substance record.
# ----------------------------------------------------------------------
class annex1record(dict):
    """
    Represents a single substance record from Annex I of EU Regulation 10/2011 EC.

    This class is a subclass of dict that stores the data for one authorized substance. It contains keys such as:
      - "record": the primary record number (1-based, as assigned from the CSV file),
      - "cid": the PubChem compound identifier (which may be None),
      - "name", "CAS", "EC", "FCM", "Ref", and additional fields such as migration limits (SML, SMLT), restrictions, and notes.

    The class provides custom __str__ and __repr__ methods for formatted display. In __repr__, only nonempty fields are shown,
    the header line begins with "record: X of Y", and text wrapping is applied so that the first line wraps at a specified width and subsequent lines wrap at a set width (with proper indenting).
    """

    def __init__(self, d, order=None, total=None):
        """
        Initialize an annex1record from dictionary d.
        order: the record number (from the CSV; if known)
        total: total number of records (if known)
        """
        if not isinstance(d,dict):
            raise TypeError("dict must be a dict not a {type(d).__name__}")
        super().__init__(d)
        # Use the value stored in the record (key "record") if available.
        self._order = d.get("record", order)
        self._total = total

    def __str__(self):
        cid = self.get("cid", None)
        order_str = f"{self._order}" if self._order is not None else "?"
        total_str = f"{self._total}" if self._total is not None else "?"
        return f"<{self.__class__.__name__} with cid:{cid} - record {order_str} of {total_str} (Annex 1 of 10/2011/EC)>"

    def __repr__(self):
        lines = []
        order_str = f"{self._order}" if self._order is not None else "?"
        total_str = f"{self._total}" if self._total is not None else "?"
        header = f" ---- [ 🇪🇺 10/2011/EC record: {order_str} of {total_str} ] ----"
        lines.append(header)
        # Define display order; note that we do not show SMLunit, SMLTunit, csvFile, or date.
        fields_order = [
            "record", "cid", "name", "CAS", "EC", "FCM", "Ref",
            "Additive_or_PPA", "Use_as_monomer_macromolecule", "FRFapplicable",
            "SML", "SMLTrestrictions", "SMLTGroupFCMsubstances", "SMLT",
            "Restriction(s)", "NotesCompliance", "Notes", "Physicalform", "ExpressedAs"
        ]
        # For each field, if value is not empty then display it.
        for key in fields_order:
            if key not in self:
                continue
            val = self[key]
            if val is None or (isinstance(val, str) and not val.strip()):
                continue
            # For SML and SMLT, append the unit inline.
            if key == "SML":
                display_val = f"{val} [mg/kg]"
            elif key == "SMLT":
                display_val = f"{val} [mg/kg]"
            else:
                display_val = str(val)
            # Wrap the value to 60 characters; subsequent lines are indented by 25 spaces.
            wrapped_val = custom_wrap(display_val, width=60, indent=" " * 22)
            lines.append(f"{key:>20}: {wrapped_val}")

        # for the extended class
        if isinstance(self,annex1record):
            lines.append(f"\n{'--- extended':>20}: properties ---\n")
            attr_order = ["cid","M","SML","SMLT","n","gFCM","gM","gMmin","gCASmin","gnamemin"]
            for attr in attr_order:
                if hasattr(self, attr):
                    val = getattr(self, attr)
                    lines.append(f"{attr:>20}: {val}")

        print("\n".join(lines))
        return str(self)

    @property
    def ispubchemok(self):
        """Returns True if the record may match in PubChem (i.e., compatible with patankar.loadpubchem)"""
        return self.get("cid") is not None and self.get("CAS") not in ("", None)

# ----------------------------------------------------------------------
# annex1record_ext: an extension of annex1 records with additional info
# The EU 10/2011 Annex 1 is not a database with unique susbtances.
# It is not fully compatible with migrant which requires a unique id.
# The extended/promoted class try to bridge them and add some info.
# We need to keep annex1=False when calling migrant to avoid circular references
# ----------------------------------------------------------------------
class annex1record_ext(annex1record):
    """extended annex1record class with additional attributes"""

    def __init__(self,rec,db=None):
        """instantiate from a record with a consolidated databse: db"""
        if not isinstance(rec,annex1record):
            raise TypeError("dict must be a annex1record not a {type(d).__name__}")
        super().__init__(rec,order=rec._order,total=rec._total)
        from patankar.loadpubchem import migrant
        if rec.ispubchemok:
            try:
                m = migrant(rec.get("CAS"),annex1=False)
                M = m.M
            except Exception:
                print(f"{rec.get('name')} could not be extended from its CAS {rec.get('CAS')}: not found")
                M = None
        else:
            M = None

        self.cid = self.get("cid")
        self.SML = self.get("SML")
        self.SMLT = self.get("SMLT")
        if self.SMLT and self.SMLT<self.SML:
            self.SML = self.SMLT
        self.gFCM = self.get("SMLTGroupFCMsubstances")
        self.n = len(self.gFCM) if self.gFCM else None
        self.M = M # molecular mass added if cid was available
        self.gM = None # molecular masses of sustances in the group
        self.gMmin = None # minimal molecular mass of substances in the group
        self.gCASmin = None # CAS of the lightest/smallest substances
        self.gnamemin = None # its name

        # we will add group info gM, gMin, gCASmin, gnamemin
        # do not forget, no bijection between rid (record id) and fcm
        # here fcm2rid[self.gFCM[0]], fcm2rid[self.gFCM[1]],...
        # are all equal, they give rids (record indices) for the group matching this FCM
        if db is not None and self.gFCM is not None:
            fcm2rid = db._fcm2rid
            # gFCM keys should be string but they can be cached as int
            if isinstance(self.gFCM[0],str):
                rids = fcm2rid[self.gFCM[0]]
            else:
                rids = fcm2rid[str(self.gFCM[0])]
            Mlist = [db._load_record(rid,db=False).M for rid in rids]
            CASlist = [db._load_record(rid,db=False).get("CAS") for rid in rids]
            namelist = [db._load_record(rid,db=False).get("name") for rid in rids]
            # Find index of the smallest non-None molecular mass
            valid_indices = [i for i, m in enumerate(Mlist) if m is not None]  # Indices of valid values
            if valid_indices:
                min_index = min(valid_indices, key=lambda i: Mlist[i])  # Index of minimum valid M value
                self.gM = [m for m in Mlist if m is not None]
                self.gMmin = Mlist[min_index]
                self.gCASmin = CASlist[min_index]
                self.gnamemin = namelist[min_index]


# ----------------------------------------------------------------------
# EuFCMannex1: the main class to manage the Annex I CSV file and cache.
# ----------------------------------------------------------------------
class EuFCMannex1:
    """
    Manages the Annex I CSV file and caches its data for efficient lookup and analysis.

    This class reads the CSV file (representing the authorized substances as provided by ECHA) and processes each row into a JSON record.
    Each record is saved in a cache directory as recXXXXX.annex1.json (where XXXXX is the 1-based record number). The class then builds a global index
    (stored as annex1_index.json in the cache folder) that maps various foreign keys—such as "name", "EC", "CAS", "FCM", "Ref", and PubChem cid—to the
    corresponding record numbers.

    Key functionality includes:
      - Retrieving records by a record number (which is now interpreted directly from the record’s "record" field, not as a Python 0-based index).
      - Lookup of records by CAS (string keys) and by PubChem cid via the bycid method.
      - Support for slicing, list-based access, callable access (e.g., dbannex1(113194) retrieves by PubChem cid or record number),
        iteration over records, and membership testing.
      - Additional search methods: byname, byEC, byCAS, byFCM, byRef, bySML, and bySMLT.
      - Handling of substances that lack a CAS number (and therefore cannot have a PubChem cid) by skipping PubChem lookup, while for substances with a nonempty CAS that fail the lookup,
        the CAS is recorded in a separate missing file ("missing.pubchem.annex1.json") to avoid future redundant queries.
      - Enhanced error messages in __getitem__ that include the invalid key and the valid range or a sample of valid keys.

    Global Rate Limiting:
      This class relies on a module-level global variable (e.g., PubChem_lastQueryTime in the patankar.loadpubchem module) to manage rate-limiting
      when performing PubChem lookups.

    Usage Example:
      >>> db = EuFCMannex1()
      >>> rec = db[1060]            # Retrieves record number 1060.
      >>> rec_by_cas = db.byCAS("102-39-6")
      >>> rec_by_cid = db.bycid(113194)
      >>> for record in db:
      ...     print(record["name"])
      >>> if 10 in db:
      ...     print("Record 10 exists.")
      >>> if 113194 in db:
      ...     print("PubChem cid 113194 exists.")
    """

    isAnnex1Initialized = False # class attribute

    def __init__(self, cache_dir="cache.EuFCMannex1", index_file="annex1_index.json", pubchem=True):
        """
        Initialize the database.

        (1) The working directory is set to the directory of __file__.
        (2) The CSV file is "fcm-and-articles-regulation--annex-i---authorised-substances-export.csv".
        (3) If the CSV file is missing, a FileNotFoundError is raised.
        (4) The class uses a cache in cache_dir and stores the global index inside that folder.
        """
        self.base_dir = os.path.dirname(__file__)
        self.csv_file = os.path.join(self.base_dir, "fcm-and-articles-regulation--annex-i---authorised-substances-export.csv")
        if not os.path.exists(self.csv_file):
            raise FileNotFoundError(f"CSV file {self.csv_file} not found.")
        self.cache_dir = os.path.join(self.base_dir, cache_dir)
        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)
        # Store the index file in the cache directory.
        self.index_file = os.path.join(self.cache_dir, index_file)
        if os.path.exists(self.index_file):
            with open(self.index_file, "r", encoding="utf-8") as f:
                self.index = json.load(f)
        else:
            self.refresh_index()
        self.order = self.index.get("order", [])
        # Internal cache for loaded records.
        self._records_cache = {}
        EuFCMannex1.isAnnex1Initialized = True

        # pubchem extensions
        # we cache the pairing fcm->rid (record id=primary key)
        self._pubchem = False # mandatory to avoid circular imports
        self._fcm2rid = self.SMLT_Groupsubstances if pubchem else None
        self._pubchem = pubchem # we enforce pubchem, the database is initialized indeed

    @classmethod
    def isindexinitialized(cls, cache_dir="cache.EuFCMannex1", index_file="annex1_index.json"):
        """Return True if the database is available"""
        return os.path.exists(os.path.join(os.path.dirname(__file__),cache_dir, index_file))

    def refresh_index(self):
        """
        Rebuild the global index by reading the CSV file and regenerating each record
        as recXXXXX.annex1.json (with record number as primary key). The index includes foreign key
        mappings for "name", "EC", "CAS", "FCM", "Ref" and a bycid index for PubChem cid.

        For substances with an empty CAS number, no PubChem lookup is performed.
        For substances with a nonempty CAS, if the lookup fails (raising ValueError) the CAS is recorded
        in a missing file ("missing.pubchem.annex1.json") so that future refreshes do not re-query PubChem.

        (Rate limiting is assumed to be managed via the module’s global variable PubChem_lastQueryTime.)
        """
        from patankar.loadpubchem import migrant # local import of migrant
        new_index = {}
        new_index["index_date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        new_index["csv_file"] = os.path.basename(self.csv_file)
        new_index["order"] = []
        for key in ["name", "EC", "CAS", "FCM", "Ref"]:
            new_index[key] = {}
        new_index["bycid"] = {}

        # Load missing CAS numbers (only for substances with a CAS)
        missing_file = os.path.join(self.cache_dir, "missing.pubchem.annex1.json")
        if os.path.exists(missing_file):
            with open(missing_file, "r", encoding="utf-8") as mf:
                missing_pubchem = json.load(mf)
        else:
            missing_pubchem = {}

        records = []
        rec_num = 0
        with open(self.csv_file, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t", quotechar='"')
            header_found = False
            for row in reader:
                if not header_found:
                    if row and row[0] == "Name":
                        header_found = True
                    continue
                if not row or len(row) < 19:
                    continue
                rec_num += 1
                rec = {}
                def yesno(val):
                    val = val.strip().lower()
                    if val == "yes":
                        return True
                    elif val == "no":
                        return False
                    return None

                # Use record number as primary key stored as "record"
                rec["record"] = rec_num
                rec["cid"] = None  # default value
                # Traceability (stored but not shown in __repr__)
                rec["csvFile"] = os.path.basename(self.csv_file)
                rec["date"] = new_index["index_date"]

                rec["name"] = row[0].strip()
                rec["EC"] = row[1].strip()
                rec["CAS"] = row[2].strip()
                rec["name2"] = row[3].strip()
                rec["CAS2"] = row[4].strip()
                try:
                    rec["FCM"] = int(row[5].strip()) if row[5].strip() not in ["", "-"] else None
                except:
                    rec["FCM"] = None
                rec["Ref"] = row[6].strip()
                rec["Additive_or_PPA"] = yesno(row[7])
                rec["Use_as_monomer_macromolecule"] = yesno(row[8])
                rec["FRFapplicable"] = yesno(row[9])

                def parse_sml(val):
                    val = val.strip()
                    if not val:
                        return SMLdefault
                    m = re.search(r"([\d\.]+)", val)
                    return float(m.group(1)) if m else SMLdefault

                rec["SML"] = parse_sml(row[10])
                rec["SMLunit"] = "mg/kg"
                rec["SMLTrestrictions"] = row[11].strip()
                val = row[12].strip()
                if val:
                    try:
                        rec["SMLTGroupFCMsubstances"] = [int(x.strip()) for x in val.split(",") if x.strip()]
                    except:
                        rec["SMLTGroupFCMsubstances"] = None
                else:
                    rec["SMLTGroupFCMsubstances"] = None
                rec["SMLT"] = parse_sml(row[13])
                if row[13].strip() == "":
                    rec["SMLT"] = rec["SML"]
                rec["SMLTunit"] = "mg/kg"
                rec["Restrictions"] = row[14].strip()
                rec["NotesCompliance"] = row[15].strip()
                rec["Notes"] = row[16].strip()
                rec["Physicalform"] = row[17].strip()
                rec["ExpressedAs"] = row[18].strip()

                # Only attempt PubChem lookup if CAS is nonempty.
                cas_val = rec["CAS"]
                if cas_val and cas_val.strip():
                    if cas_val in missing_pubchem:
                        cid_val = missing_pubchem[cas_val]
                    else:
                        try:
                            cid_val = migrant(cas_val,annex1=False).cid
                        except ValueError:
                            printWARN(f"🇪🇺 Warning: substance {rec['name']} (CAS {cas_val}) not found in PubChem.")
                            cid_val = None
                            missing_pubchem[cas_val] = None
                else:
                    # No CAS provided, so we know there will be no cid.
                    cid_val = None
                rec["cid"] = cid_val

                # Order the record keys: first "record", then "cid", then others.
                ordered_rec = {
                    "record": rec["record"],
                    "cid": rec["cid"],
                    "name": rec["name"],
                    "CAS": rec["CAS"],
                    "EC": rec["EC"],
                    "FCM": rec["FCM"],
                    "Ref": rec["Ref"],
                    "Additive_or_PPA": rec["Additive_or_PPA"],
                    "Use_as_monomer_macromolecule": rec["Use_as_monomer_macromolecule"],
                    "FRFapplicable": rec["FRFapplicable"],
                    "SML": rec["SML"],
                    "SMLunit": rec["SMLunit"],
                    "SMLTrestrictions": rec["SMLTrestrictions"],
                    "SMLTGroupFCMsubstances": rec["SMLTGroupFCMsubstances"],
                    "SMLT": rec["SMLT"],
                    "SMLTunit": rec["SMLTunit"],
                    "Restrictions": rec["Restrictions"],
                    "NotesCompliance": rec["NotesCompliance"],
                    "Notes": rec["Notes"],
                    "Physicalform": rec["Physicalform"],
                    "ExpressedAs": rec["ExpressedAs"],
                    "engine": f"SFPPy: {os.path.basename(__file__)}",
                    "csvFile": rec["csvFile"],
                    "date": rec["date"]
                }

                rec_filename = f"rec{rec_num:05d}.annex1.json"
                json_filename = os.path.join(self.cache_dir, rec_filename)
                with open(json_filename, "w", encoding="utf-8") as jf:
                    json.dump(ordered_rec, jf, ensure_ascii=False, indent=2)

                new_index["order"].append(rec_num)
                for key in ["name", "EC", "CAS", "FCM", "Ref"]:
                    value = rec[key]
                    if value not in new_index[key]:
                        new_index[key][value] = []
                    new_index[key][value].append(rec_num)
                if cid_val is not None:
                    new_index["bycid"][cid_val] = rec_num
                records.append(ordered_rec)
        with open(self.index_file, "w", encoding="utf-8") as f:
            json.dump(new_index, f, ensure_ascii=False, indent=2)
        with open(missing_file, "w", encoding="utf-8") as mf:
            json.dump(missing_pubchem, mf, ensure_ascii=False, indent=2)
        self.index = new_index
        self.order = new_index.get("order", [])
        self._records_cache = {}

    def _load_record(self, rec_id, order=None, db=False):
        """
        Load a record (as an annex1record) from its cached JSON file.
        If the file does not exist, return None with a warning.
        Extended records are managed via the global flag self._pubchem
        The local flag db sets whether the record will be informed or not from the full database
        """
        if rec_id in self._records_cache:
            if self._pubchem:
                if db:
                    return annex1record_ext(self._records_cache[rec_id],self)
                else:
                    return annex1record_ext(self._records_cache[rec_id])
            else:
                return self._records_cache[rec_id]
        json_filename = os.path.join(self.cache_dir, f"rec{rec_id:05d}.annex1.json")
        if not os.path.exists(json_filename):
            printWARN(f"🇪🇺 Warning: Record file for record {rec_id} not found.")
            return None
        with open(json_filename, "r", encoding="utf-8") as jf:
            rec = json.load(jf)
        # Pass the record's own "record" field as order.
        record_obj = annex1record(rec, order=rec.get("record"), total=len(self.order))
        self._records_cache[rec_id] = record_obj
        if self._pubchem:
            self._pubchem = False # mandatory to avoid circular imports
            self._fcm2rid = self.SMLT_Groupsubstances # we refresh the cache with the new substance
            self._pubchem = True
            if db:
                return annex1record_ext(record_obj,self)
            else:
                return annex1record_ext(record_obj)
        else:
            return record_obj

    def __getitem__(self, key):
        """
        __getitem__ supports:
         - integer keys: interpreted as record numbers.
           If the record number is not found, an error is raised showing the entered value and the valid range.
         - slices: returns a list of records whose record numbers fall within the slice.
         - list/tuple: returns a list of corresponding records.
         - string keys: interpreted as CAS numbers.
           If the CAS number is not found, an error is raised including the entered value and a sample of valid keys.
        """
        if isinstance(key, slice):
            # Interpret the slice start and stop as record numbers.
            start = key.start if key.start is not None else min(self.order)
            stop = key.stop if key.stop is not None else max(self.order) + 1
            rec_ids = [rid for rid in self.order if start <= rid < stop]
            if not rec_ids:
                raise KeyError(f"No records found in range {start} to {stop - 1}. Valid record numbers range from {min(self.order)} to {max(self.order)}.")
            return [self._load_record(rid, order=rid) for rid in rec_ids]
        elif isinstance(key, int):
            # Interpret the integer key directly as a record number.
            if key in self.order:
                return self._load_record(key, order=key)
            else:
                raise KeyError(f"Record number {key} not found. Valid record numbers range from {min(self.order)} to {max(self.order)}.")
        elif isinstance(key, (list, tuple)):
            return [self.__getitem__(k) for k in key]
        elif isinstance(key, str):
            # Interpret a string key as a CAS number.
            if key in self.index.get("CAS", {}):
                rec_ids = self.index["CAS"][key]
                if len(rec_ids) == 1:
                    return self._load_record(rec_ids[0], order=rec_ids[0])
                else:
                    return [self._load_record(rid, order=rid) for rid in rec_ids]
            else:
                available = list(self.index.get("CAS", {}).keys())
                sample = ", ".join(available[:10]) + (" ..." if len(available) > 10 else "")
                raise KeyError(f"CAS key '{key}' not found in index. Valid CAS keys include: {sample}")
        else:
            raise KeyError(f"Unsupported key type: {type(key)}")


    def __call__(self, *args):
        """
        Callable access. For example:
         - dbannex1(cid) returns the record for a given PubChem cid (via bycid),
         - dbannex1(rec) returns the record for a given record number,
         - dbannex1(cid1, cid2, ...) or dbannex1([cid1, cid2, ...]) returns a list.
         Strings are interpreted as CAS numbers.
        """
        if len(args) == 1 and isinstance(args[0], (list, tuple)):
            args = args[0]
        results = []
        for arg in args:
            if isinstance(arg, int):
                argkey = str(arg)
                if arg in self.order:
                    results.append(self._load_record(arg))
                elif "bycid" in self.index and argkey in self.index["bycid"]:
                    rec_id = self.index["bycid"][argkey]
                    results.append(self._load_record(rec_id))
                else:
                    printWARN(f"🇪🇺 Warning: Record for identifier {arg} not found.")
                    results.append(None)
            elif isinstance(arg, str):
                result_item = self.__getitem__(arg)
                if isinstance(result_item, list):
                    results.extend(result_item)
                else:
                    results.append(result_item)
            else:
                raise KeyError(f"Unsupported key type in call: {type(arg)}")
        if len(results) == 1:
            return results[0]
        return results

    # --------------------------
    # Additional search methods
    # --------------------------
    def byname(self, name):
        name = name[0] if isinstance(name,list) else name
        rec_ids = self.index.get("name", {}).get(name, [])
        return [self._load_record(rid, order=rid) for rid in rec_ids]

    def byEC(self, ec):
        ec = ec[0] if isinstance(ec,list) else ec
        rec_ids = self.index.get("EC", {}).get(ec, [])
        return [self._load_record(rid, order=rid) for rid in rec_ids]

    def byCAS(self, cas):
        cas = cas[0] if isinstance(cas,list) else cas
        rec_ids = self.index.get("CAS", {}).get(cas, [])
        if len(rec_ids) == 1:
            return self._load_record(rec_ids[0], order=rec_ids[0],db=True)
        else:
            return [self._load_record(rid, order=rid,db=True) for rid in rec_ids]

    def byFCM(self, fcm):
        fcm = fcm[0] if isinstance(fcm,list) else fcm
        fcm = str(fcm) if isinstance(fcm,int) else fcm
        rec_ids = self.index.get("FCM", {}).get(fcm, [])
        return [self._load_record(rid, order=rid) for rid in rec_ids]

    def byRef(self, ref):
        ref = ref[0] if isinstance(ref,list) else ref
        rec_ids = self.index.get("Ref", {}).get(ref, [])
        return [self._load_record(rid, order=rid) for rid in rec_ids]

    def bycid(self, cid, verbose=True):
        """
        Search for a record by PubChem cid.
        """
        cid = cid[0] if isinstance(cid,list) else cid
        cidkey = str(cid)
        if "bycid" in self.index and cidkey in self.index["bycid"]:
            rec_id = self.index["bycid"][str(cid)]
            return self._load_record(rec_id, order=rec_id,db=True)
        else:
            if verbose:
                printWARN(f"⚠️ Warning: No 🇪🇺 10/2011/EC record found for PubChem cid {cid}.")
            return None

    def bySML(self, min_val, max_val):
        results = []
        for rid in self.order:
            rec = self._load_record(rid)
            if rec is not None:
                sml = rec.get("SML", SMLdefault)
                if min_val <= sml <= max_val:
                    results.append(rec)
        return results

    def bySMLT(self, min_val, max_val):
        results = []
        for rid in self.order:
            rec = self._load_record(rid)
            if rec is not None:
                smlt = rec.get("SMLT", SMLdefault)
                if min_val <= smlt <= max_val:
                    results.append(rec)
        return results

    def __iter__(self):
        for rid in self.order:
            yield self._load_record(rid, order=rid)

    def __len__(self):
        return len(self.order)

    def __contains__(self, item):
        if isinstance(item,list):
            item = item[0]
        if isinstance(item, int):
            argkey = str(item)
            return item in self.order or ("bycid" in self.index and argkey in self.index["bycid"])
        if isinstance(item, str):
            return item in self.index.get("CAS", {})
        return False

    @property
    def SMLT_Groupsubstances(self):
        """
            Returns a dictionary G so that G[fcm] = [rid1, rid2] (cached in memory and on disk)

            - **Hidden Field Cache:**
            The property first checks if `self._SMLT_Groupsubstances` exists and returns it immediately if so.

            - **File Cache:**
            The cache file path is built by taking `self.index_file`, splitting it into a base and extension, then appending “.group” before the extension.
            If this cache file exists, the code attempts to load it. Any issues during file reading (e.g., corruption) will fall back to recomputing the groups.

            - **Recomputation and Write-back:**
            If the cache file does not exist or fails to load, the code computes the groups dictionary, caches it in memory, and writes the result to the cache file.

            This double caching strategy ensures that once the groups are computed, subsequent accesses are fast (both from memory and from disk), while still allowing a persistent cache that survives between sessions.
        """
        # If already cached in the hidden field, return it immediately.
        if hasattr(self, '_SMLT_Groupsubstances'):
            return self._SMLT_Groupsubstances

        # Build the cache file name using self.index_file:
        base, ext = os.path.splitext(self.index_file)
        cache_file = base + ".group" + ext  # e.g., annex1_index.group.json

        # Try to load the cache from the file if it exists.
        if os.path.exists(cache_file):
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    groups = json.load(f)
                self._SMLT_Groupsubstances = groups
                return groups
            except Exception as e:
                printWARN(f"🇪🇺 Warning: Failed to load group cache from {cache_file}: {e}")

        # Compute the groups if not cached.
        groups = {}
        for rid in self.order:
            rec = self._load_record(rid)
            lst = rec.get("SMLTGroupFCMsubstances")
            if lst:
                for fcm in lst:
                    groups.setdefault(fcm, []).append(rec.get("record"))

        # Cache the result in the hidden field.
        self._SMLT_Groupsubstances = groups

        # Write the computed groups to the cache file.
        try:
            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(groups, f, indent=2, ensure_ascii=False)
        except Exception as e:
            printWARN(f"🇪🇺 Warning: Failed to write group cache to {cache_file}: {e}")

        return groups


    def __repr__(self):
        csv_filename = os.path.basename(self.csv_file)
        index_date = self.index.get("index_date", "unknown")
        print(f"Annex 1 of 🇪🇺EU regulation 10/2011 EC ({len(self.order)} records)")
        print(f"Imported from CSV {csv_filename} and indexed on {index_date}")
        return str(self)

    def __str__(self):
        return f"<{self.__class__.__name__}: {len(self.order)} records (Annex 1 of 10/2011/EC)>"

# -------------------------------------------------------------------
# Example usage (for debugging / standalone tests)
# -------------------------------------------------------------------
if __name__ == "__main__":
    dbannex1 = EuFCMannex1(pubchem=True) # we promote the whole database
    rec = dbannex1.bycid(6581)
    repr(rec)

    print(repr(dbannex1))
    print(str(dbannex1))
    first_record = dbannex1[1]
    print(first_record)

    # Search methods examples:
    rec_by_cas = dbannex1.byCAS("102-39-6")
    rec_by_name = dbannex1.byname("(1,3-phenylenedioxy)diacetic acid")

    # Search by PubChem cid:
    rec_by_cid = dbannex1.bycid(113194)  # if valid
    rec_by_cid_missing = dbannex1.bycid(1023456)  # will warn and return None

    # Search by SML range:
    sml_records = dbannex1.bySML(0.01, 1.0)

    # Callable access:
    rec = dbannex1(113194)         # by cid (if exists) or by record number
    rec_list = dbannex1(1023456, "102-39-6")

    # Iterate over records:
    for rec in dbannex1:
        print(rec.get("cid"))

    # Test membership:
    if 10 in dbannex1:
        print("Record 10 exists.")
    if 113194 in dbannex1:
        print("PubChem cid 113194 exists.")
