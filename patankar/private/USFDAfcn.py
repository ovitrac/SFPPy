
__all__ = ['USFDAfcn', 'clean_html', 'custom_wrap', 'fcnrecord', 'fcnrecord_ext', 'parse_food_contact_substance', 'printWARN']
"""
Module: USFDAfcn

This module implements a robust database manager for the US Inventory of Effective Food Contact
Substance (FCS) notifications as provided by the US FDA. The CSV file (exported from the US FDA website)
contains the official Food Contact Substance Notification list but with information not well organized.

Each row in the CSV file is converted into a JSON record (named as recXXXXX.fcn.json in the cache directory)
with the following ordered fields:
  - "record": the sequential record number (1-based)
  - "cid": the PubChem compound identifier (or list of cids if a mixture)
  - "name": the chemical name (a string for single substances or a list for mixtures)
  - "CAS": the CAS number (a string or list if a mixture)
  - "FCNNo": the FCN number extracted from the first column (e.g. =T("2355") becomes "2355")
  - "FoodContactSubstance": the full original field text from the CSV
  - "mixture": boolean (True if the substance is a mixture, i.e. compound and CAS are lists)
  - "FCSreplacedby": field as provided (if the notification has been replaced)
  - "FCSreplacedby_record": the record number corresponding to the replacement (if available)
  - "notifier": the notifier field from the CSV
  - "manufacturer": from the CSV field "Manufacture/Supplier"
  - "NotificationDate": the notification date from the CSV

At the end, traceability fields "engine", "csfile", and "date" are added.

The module also builds a global index (stored in fcn_index.json) mapping key fields:
  - "name", "CAS", "FCNNo", "FCSreplacedby" and PubChem cid ("bycid")
to the corresponding record numbers.

Key features include:
  - Robust parsing of the FoodContactSubstance field using regular expressions.
  - Handling of mixtures (multiple chemicals in one notification) where the chemical names and CAS numbers
    are stored as lists.
  - PubChem lookup for each CAS number (implemented via patankar.loadpubchem.migrant) to retrieve the corresponding cid.
  - A secondary pass to resolve FCSreplacedby: if a record indicates that it has been replaced by another notification,
    a new field "FCSreplacedby_record" is added corresponding to the record number of the replacement notification.
  - Lookup methods by record number, by CAS, by FCNNo, and by PubChem cid.
  - Caching of missing PubChem results in "missing.pubchem.fcn.json" to avoid repeated failed queries.

@version: 1.41
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2024-01-10
@rev: 2025-03-31

"""

import os, csv, json, datetime, time, re, textwrap

__project__ = "SFPPy"
__author__ = "Olivier Vitrac"
__copyright__ = "Copyright 2022"
__credits__ = ["Olivier Vitrac"]
__license__ = "MIT"
__maintainer__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@agroparistech.fr"
__version__ = "1.41"


# Default PubChem lookup error value (if not found, value remains None)
# (This can be adapted if needed)
DEFAULT_PUBCHEM = None

# Module-level variables to track last warning message and its timestamp
_LAST_WARN_ = None
_T_LAST_WARN_ = 0.0

# ----------------------------------------------------------------------
# Utility function: custom_wrap (for pretty printing)
# ----------------------------------------------------------------------
def custom_wrap(text, width=60, indent=" " * 22):
    # Wrap first line to specified width.
    first_line = textwrap.wrap(text, width=width)
    if not first_line:
        return ""
    first = first_line[0]
    remaining = text[len(first):].lstrip()
    subsequent_lines = textwrap.wrap(remaining, width=width)
    wrapped = [first] + [indent + line for line in subsequent_lines]
    return "\n".join(wrapped)

# ----------------------------------------------------------------------
# Show warnings without repeating them
# ----------------------------------------------------------------------
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

# ----------------------------------------------------------------------
# *New helper function to remove HTML tags*
# ----------------------------------------------------------------------
def clean_html(text):
    return re.sub(r'<[^>]*>', '', text).strip()

# ----------------------------------------------------------------------
# Helper function to parse FoodContactSubstance field
# ----------------------------------------------------------------------
def parse_food_contact_substance(text):
    """
    Parses the full FoodContactSubstance field to extract the chemical name(s) and CAS number(s).

    For a single chemical, the name is defined as the portion of the text from the start until the first
    occurrence of "(CAS Reg. No." (after removing a trailing phrase such as " (produced by ...").

    For mixtures (multiple occurrences), a fallback using re.findall is used.

    Returns:
       name: a string if one match is found, or a list of names for mixtures.
       CAS: a string if one match is found, or a list of CAS numbers for mixtures.
       mixture: Boolean, True if more than one match is found.
    """
    text = text.strip()
    # If there is no CAS pattern, return the full text.
    if "(CAS Reg. No." not in text:
        return text, "", False

    count = text.count("(CAS Reg. No.")
    if count == 1:
        idx = text.find("(CAS Reg. No.")
        candidate = text[:idx].strip()
        # If the candidate contains a trailing phrase like " (produced by", remove it.
        p_idx = candidate.lower().rfind(" (produced by")
        if p_idx != -1:
            candidate = candidate[:p_idx].strip()
        candidate = candidate.rstrip(",;").strip()
        m = re.search(r'\(CAS Reg\. No\. ([\d-]+)\)', text[idx:])
        cas = m.group(1) if m else ""
        return candidate, cas.strip(), False
    else:
        # For mixtures, use re.findall to capture all occurrences.
        # This pattern attempts to capture a chemical name (without leading commas/semicolons)
        # and its corresponding CAS number.
        pattern = r'([^,(;]+(?:\([^)]*\))?[^,(;]*)\s*\(CAS Reg\. No\. ([\d-]+)\)'
        matches = re.findall(pattern, text)
        if not matches:
            # Fallback: treat as single if no matches are found.
            idx = text.find("(CAS Reg. No.")
            candidate = text[:idx].strip().rstrip(",;")
            m = re.search(r'\(CAS Reg\. No\. ([\d-]+)\)', text[idx:])
            cas = m.group(1) if m else ""
            return candidate, cas.strip(), False
        names = []
        cas_list = []
        for nm, cas in matches:
            nm_clean = nm.strip().lstrip(",;").strip()
            names.append(nm_clean)
            cas_list.append(cas.strip())
        if len(names) == 1:
            return names[0], cas_list[0], False
        return names, cas_list, True

# ----------------------------------------------------------------------
# Class: fcnrecord
# ----------------------------------------------------------------------
class fcnrecord(dict):
    """
    Represents a single Food Contact Substance Notification record from the US FDA database.

    Keys include:
      - "record": the sequential record number (1-based)
      - "cid": PubChem compound identifier (or list of cids for mixtures; may be None)
      - "name": interpreted chemical name (string or list for mixtures)
      - "CAS": CAS number (string or list for mixtures)
      - "FCNNo": the FCN number extracted from the CSV row (e.g. from =T("2355"))
      - "FoodContactSubstance": the full field text from the CSV
      - "mixture": boolean flag indicating if the record represents a mixture
      - "FCSreplacedby": original field (as provided) indicating a replacement notification
      - "FCSreplacedby_record": record number (if found) corresponding to the replacement notification
      - "notifier", "manufacturer", "NotificationDate": additional fields from the CSV
      - Traceability fields: "engine", "csfile", "date"
    """

    def __init__(self, d, order=None, total=None):
        if not isinstance(d, dict):
            raise TypeError("Input must be a dict, not a {}".format(type(d).__name__))
        super().__init__(d)
        self._order = d.get("record", order)
        self._total = total

    def __str__(self):
        cid = self.get("cid", None)
        order_str = f"{self._order}" if self._order is not None else "?"
        total_str = f"{self._total}" if self._total is not None else "?"
        return f"<{self.__class__.__name__} with cid:{cid} - record {order_str} of {total_str} (US FDA FCS)>"

    def __repr__(self):
        lines = []
        order_str = f"{self._order}" if self._order is not None else "?"
        total_str = f"{self._total}" if self._total is not None else "?"
        header = f" ---- [ US FDA FCS record: {order_str} of {total_str} ] ----"
        lines.append(header)
        fields_order = [
            "record", "cid", "name", "CAS", "FCNNo",
            "FoodContactSubstance", "mixture", "FCSreplacedby", "FCSreplacedby_record",
            "notifier", "manufacturer", "NotificationDate"
        ]
        for key in fields_order:
            if key not in self:
                continue
            val = self[key]
            if val is None or (isinstance(val, str) and not val.strip()):
                continue
            wrapped_val = custom_wrap(str(val), width=60, indent=" " * 22)
            lines.append(f"{key:>20}: {wrapped_val}")
        for key in ["engine", "csfile", "date"]:
            if key in self:
                wrapped_val = custom_wrap(str(self[key]), width=60, indent=" " * 22)
                lines.append(f"{key:>20}: {wrapped_val}")
        return "\n".join(lines)

    @property
    def ispubchemok(self):
        cas = self.get("CAS")
        if self.get("mixture", False):
            return bool(cas and any(c.strip() for c in cas))
        return cas not in ("", None)

# ----------------------------------------------------------------------
# Class: fcnrecord_ext
# ----------------------------------------------------------------------
class fcnrecord_ext(fcnrecord):
    """
    Extended fcnrecord that automatically retrieves additional chemical information from PubChem.

    For each CAS number (or each CAS in a mixture) the PubChem lookup is performed.
    The field "cid" is updated to be either a single PubChem CID or a list of CIDs.
    """
    def __init__(self, rec, db=None, verbosity=False):
        """
        Instantiate from a base fcnrecord.
        If a valid CAS is available, perform PubChem lookup via the 'migrant' function.
        """
        if not isinstance(rec, fcnrecord):
            raise TypeError("Input must be an fcnrecord, not a {}".format(type(rec).__name__))
        super().__init__(rec, order=rec._order, total=rec._total)
        from patankar.loadpubchem import migrant
        if self.ispubchemok:
            if self.get("mixture", False):
                cids = []
                for cas in self.get("CAS", []):
                    try:
                        m = migrant(cas, annex1=False)
                        cids.append(m.cid)
                    except Exception:
                        if verbosity:
                            printWARN(f"🇺🇸 Warning: PubChem lookup failed for CAS {cas} in record {self.get('record')}")
                        cids.append(DEFAULT_PUBCHEM)
                self.cid = cids
            else:
                cas = self.get("CAS")
                try:
                    m = migrant(cas, annex1=False)
                    self.cid = m.cid
                except Exception:
                    if verbosity:
                        printWARN(f"🇺🇸 Warning: PubChem lookup failed for CAS {cas} in record {self.get('record')}")
                    self.cid = DEFAULT_PUBCHEM
        else:
            self.cid = None


# ----------------------------------------------------------------------
# Class: USFDAfcn
# ----------------------------------------------------------------------
class USFDAfcn:
    """
    Manages the US FDA Food Contact Substance Notification CSV file and caches its data for efficient lookup.

    This class reads the official CSV file and processes each row into a JSON record (recXXXXX.fcn.json)
    stored in a cache directory (default "cache.USFDAfcn"). It then builds a global index (stored as fcn_index.json)
    mapping key fields—such as "name", "CAS", "FCNNo", "FCSreplacedby" and PubChem cid—to their corresponding record numbers.

    Lookup methods include:
      - __getitem__: lookup by sequential record number or by CAS string.
      - __call__: callable access by record number, PubChem cid, or CAS.
      - byname, byCAS, byFCNNo, byFCSreplacedby, bycid: dedicated search methods.

    Global Rate Limiting:
      PubChem lookups are assumed to be managed by a module-level variable in patankar.loadpubchem.
    """
    isInitialized = False

    def __init__(self, cache_dir="cache.USFDAfcn", index_file="fcn_index.json", pubchem=True):
        self.base_dir = os.path.dirname(__file__)
        self.csv_file = os.path.join(self.base_dir, "FCN.csv")
        if not os.path.exists(self.csv_file):
            raise FileNotFoundError(f"CSV file {self.csv_file} not found.")
        self.cache_dir = os.path.join(self.base_dir, cache_dir)
        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)
        self.index_file = os.path.join(self.cache_dir, index_file)
        if os.path.exists(self.index_file):
            with open(self.index_file, "r", encoding="utf-8") as f:
                self.index = json.load(f)
        else:
            self.refresh_index()
        self.order = self.index.get("order", [])
        self._records_cache = {}
        self._pubchem = pubchem
        USFDAfcn.isInitialized = True

    @classmethod
    def isindexinitialized(cls, cache_dir="cache.USFDAfcn", index_file="fcn_index.json"):
        return os.path.exists(os.path.join(os.path.dirname(__file__), cache_dir, index_file))

    def refresh_index(self):
        from patankar.loadpubchem import migrant
        new_index = {}
        index_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        new_index["index_date"] = index_date
        new_index["csv_file"] = os.path.basename(self.csv_file)
        new_index["order"] = []
        for key in ["name", "CAS", "FCNNo", "FCSreplacedby"]:
            new_index[key] = {}
        new_index["bycid"] = {}

        missing_file = os.path.join(self.cache_dir, "missing.pubchem.fcn.json")
        if os.path.exists(missing_file):
            with open(missing_file, "r", encoding="utf-8") as mf:
                missing_pubchem = json.load(mf)
        else:
            missing_pubchem = {}

        records_list = []
        fcnno2recid = {}

        rec_num = 0
        with open(self.csv_file, "r", encoding="latin1") as f:
            while True:
                pos = f.tell()
                line = f.readline()
                if not line:
                    break
                if line.lstrip().startswith("FCN No"):
                    f.seek(pos)
                    break
            reader = csv.reader(f, delimiter=",")
            header = next(reader, None)
            for row in reader:
                if not row or len(row) < 7:
                    continue
                rec_num += 1
                rec = {}
                fcn_no_match = re.search(r'=T\("(\d+)"\)', row[0])
                if fcn_no_match:
                    fcn_no = fcn_no_match.group(1)
                else:
                    fcn_no = row[0].strip()
                rec["FCNNo"] = fcn_no

                # Parse Food Contact Substance from column 6
                fcs_text = row[6].strip()
                rec["FoodContactSubstance"] = fcs_text
                name_val, cas_val, is_mixture = parse_food_contact_substance(fcs_text)
                rec["name"] = name_val
                rec["CAS"] = cas_val
                rec["mixture"] = is_mixture

                # FCS REPLACED BY from column 7
                rec["FCSreplacedby"] = row[7].strip() if len(row) > 7 else ""
                # Notifier from column 8
                rec["notifier"] = row[8].strip() if len(row) > 8 else ""
                # *Clean manufacturer field by removing HTML tags (from column 9)*
                rec["manufacturer"] = clean_html(row[9].strip()) if len(row) > 9 else ""
                # Notification Date from column 12 (Effective Date)
                rec["NotificationDate"] = row[12].strip() if len(row) > 12 else ""

                rec["csvFile"] = os.path.basename(self.csv_file)
                rec["date"] = index_date
                rec["engine"] = f"SFPPy: {os.path.basename(__file__)}"
                rec["record"] = rec_num
                rec["cid"] = None

                if rec["CAS"] and not rec.get("mixture", False):
                    cas = rec["CAS"]
                    if cas in missing_pubchem:
                        cid_val = missing_pubchem[cas]
                    else:
                        try:
                            cid_val = migrant(cas, annex1=False).cid
                        except Exception:
                            printWARN(f"🇺🇸 Warning: PubChem lookup failed for {rec['name']} (CAS {cas}).")
                            cid_val = None
                            missing_pubchem[cas] = None
                    rec["cid"] = cid_val
                elif rec["CAS"] and rec.get("mixture", False):
                    cid_list = []
                    for cas in rec["CAS"]:
                        if cas in missing_pubchem:
                            cid_val = missing_pubchem[cas]
                        else:
                            try:
                                cid_val = migrant(cas, annex1=False).cid
                            except Exception:
                                printWARN(f"🇺🇸 Warning: PubChem lookup failed for component with CAS {cas} in record {rec_num}.")
                                cid_val = None
                                missing_pubchem[cas] = None
                        cid_list.append(cid_val)
                    rec["cid"] = cid_list

                ordered_rec = {
                    "record": rec["record"],
                    "cid": rec["cid"],
                    "name": rec["name"],
                    "CAS": rec["CAS"],
                    "FCNNo": rec["FCNNo"],
                    "FoodContactSubstance": rec["FoodContactSubstance"],
                    "mixture": rec["mixture"],
                    "FCSreplacedby": rec["FCSreplacedby"],
                    "FCSreplacedby_record": None,
                    "notifier": rec["notifier"],
                    "manufacturer": rec["manufacturer"],
                    "NotificationDate": rec["NotificationDate"],
                    "engine": rec["engine"],
                    "csfile": rec["csvFile"],
                    "date": rec["date"]
                }
                rec_filename = f"rec{rec_num:05d}.fcn.json"
                json_filename = os.path.join(self.cache_dir, rec_filename)
                with open(json_filename, "w", encoding="utf-8") as jf:
                    json.dump(ordered_rec, jf, ensure_ascii=False, indent=2)

                new_index["order"].append(rec_num)
                if rec.get("mixture", False):
                    for nm in rec["name"]:
                        new_index["name"].setdefault(nm, []).append(rec_num)
                else:
                    new_index["name"].setdefault(rec["name"], []).append(rec_num)
                if rec.get("mixture", False):
                    for cas in rec["CAS"]:
                        new_index["CAS"].setdefault(cas, []).append(rec_num)
                else:
                    new_index["CAS"].setdefault(rec["CAS"], []).append(rec_num)
                new_index["FCNNo"].setdefault(rec["FCNNo"], []).append(rec_num)
                if rec["FCSreplacedby"]:
                    new_index["FCSreplacedby"].setdefault(rec["FCSreplacedby"], []).append(rec_num)
                if rec["cid"] is not None:
                    if rec.get("mixture", False):
                        for cid in rec["cid"]:
                            if cid is not None:
                                new_index["bycid"][str(cid)] = rec_num
                    else:
                        new_index["bycid"][str(rec["cid"])] = rec_num

                fcnno2recid[rec["FCNNo"]] = rec_num
                records_list.append((rec_num, ordered_rec))

        # Second pass: resolve FCSreplacedby_record using FCNNo mapping
        for rec_id, rec in records_list:
            fcsrep_field = rec.get("FCSreplacedby", "").strip()
            if fcsrep_field:
                # *Extract FCN number from FCSreplacedby field if present*
                fcsrep_match = re.search(r'FCN\s*(\d+)', fcsrep_field, re.IGNORECASE)
                if fcsrep_match:
                    fcsrep = fcsrep_match.group(1)
                else:
                    fcsrep = fcsrep_field
                rep_recid = fcnno2recid.get(fcsrep)
                rec["FCSreplacedby_record"] = rep_recid
                json_filename = os.path.join(self.cache_dir, f"rec{rec_id:05d}.fcn.json")
                try:
                    with open(json_filename, "w", encoding="utf-8") as jf:
                        json.dump(rec, jf, ensure_ascii=False, indent=2)
                except Exception as e:
                    printWARN(f"🇺🇸 Warning: Could not update FCSreplacedby_record in {json_filename}: {e}")

        with open(self.index_file, "w", encoding="utf-8") as f_index:
            json.dump(new_index, f_index, ensure_ascii=False, indent=2)
        with open(missing_file, "w", encoding="utf-8") as mf:
            json.dump(missing_pubchem, mf, ensure_ascii=False, indent=2)
        self.index = new_index
        self.order = new_index.get("order", [])
        self._records_cache = {}

    def _load_record(self, rec_id, order=None, db=False):
        if rec_id in self._records_cache:
            if self._pubchem:
                if db:
                    return fcnrecord_ext(self._records_cache[rec_id], self)
                else:
                    return fcnrecord_ext(self._records_cache[rec_id])
            else:
                return self._records_cache[rec_id]
        json_filename = os.path.join(self.cache_dir, f"rec{rec_id:05d}.fcn.json")
        if not os.path.exists(json_filename):
            printWARN(f"🇺🇸 Warning: Record file for record {rec_id} not found.")
            return None
        with open(json_filename, "r", encoding="utf-8") as jf:
            rec = json.load(jf)
        record_obj = fcnrecord(rec, order=rec.get("record"), total=len(self.order))
        self._records_cache[rec_id] = record_obj
        if self._pubchem:
            return fcnrecord_ext(record_obj, self)
        else:
            return record_obj

    def __getitem__(self, key):
        if isinstance(key, slice):
            start = key.start if key.start is not None else min(self.order)
            stop = key.stop if key.stop is not None else max(self.order) + 1
            rec_ids = [rid for rid in self.order if start <= rid < stop]
            if not rec_ids:
                raise KeyError(f"No records found in range {start} to {stop - 1}. Valid records range from {min(self.order)} to {max(self.order)}.")
            return [self._load_record(rid, order=rid) for rid in rec_ids]
        elif isinstance(key, int):
            if key in self.order:
                return self._load_record(key, order=key)
            else:
                raise KeyError(f"Record number {key} not found. Valid records range from {min(self.order)} to {max(self.order)}.")
        elif isinstance(key, (list, tuple)):
            return [self.__getitem__(k) for k in key]
        elif isinstance(key, str):
            # First, try CAS lookup.
            if key in self.index.get("CAS", {}):
                rec_ids = self.index["CAS"][key]
                if len(rec_ids) == 1:
                    return self._load_record(rec_ids[0], order=rec_ids[0])
                else:
                    return [self._load_record(rid, order=rid) for rid in rec_ids]
            # Then, if key is all digits, try FCNNo lookup.
            elif key.isdigit() and key in self.index.get("FCNNo", {}):
                rec_ids = self.index["FCNNo"][key]
                if len(rec_ids) == 1:
                    return self._load_record(rec_ids[0], order=rec_ids[0])
                else:
                    return [self._load_record(rid, order=rid) for rid in rec_ids]
            else:
                available = list(self.index.get("CAS", {}).keys()) + list(self.index.get("FCNNo", {}).keys())
                sample = ", ".join(available[:10]) + (" ..." if len(available) > 10 else "")
                raise KeyError(f"Key '{key}' not found in index. Valid keys include: {sample}")
        else:
            raise KeyError(f"Unsupported key type: {type(key)}")

    def __call__(self, *args):
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
                    printWARN(f"🇺🇸 Warning: Record for identifier {arg} not found.")
                    results.append(None)
            elif isinstance(arg, str):
                result_item = self.__getitem__(arg)
                if isinstance(result_item, list):
                    results.extend(result_item)
                else:
                    results.append(result_item)
            else:
                raise KeyError(f"Unsupported key type in call: {type(arg)}")
        return results[0] if len(results) == 1 else results

    def byname(self, name):
        name = name[0] if isinstance(name, list) else name
        rec_ids = self.index.get("name", {}).get(name, [])
        return [self._load_record(rid, order=rid) for rid in rec_ids]

    def byCAS(self, cas):
        cas = cas[0] if isinstance(cas, list) else cas
        rec_ids = self.index.get("CAS", {}).get(cas, [])
        if len(rec_ids) == 1:
            return self._load_record(rec_ids[0], order=rec_ids[0])
        else:
            return [self._load_record(rid, order=rid) for rid in rec_ids]

    def byFCNNo(self, fcn_no):
        fcn_no = fcn_no[0] if isinstance(fcn_no, list) else fcn_no
        rec_ids = self.index.get("FCNNo", {}).get(fcn_no, [])
        if len(rec_ids) == 1:
            return self._load_record(rec_ids[0], order=rec_ids[0])
        else:
            return [self._load_record(rid, order=rid) for rid in rec_ids]

    def byFCSreplacedby(self, fcsrep):
        fcsrep = fcsrep[0] if isinstance(fcsrep, list) else fcsrep
        rec_ids = self.index.get("FCSreplacedby", {}).get(fcsrep, [])
        if len(rec_ids) == 1:
            return self._load_record(rec_ids[0], order=rec_ids[0])
        else:
            return [self._load_record(rid, order=rid) for rid in rec_ids]

    def bycid(self, cid, verbose=True):
        cid = cid[0] if isinstance(cid, list) else cid
        cidkey = str(cid)
        if "bycid" in self.index and cidkey in self.index["bycid"]:
            rec_id = self.index["bycid"][cidkey]
            return self._load_record(rec_id, order=rec_id)
        else:
            if verbose:
                printWARN(f"⚠️ Warning: No 🇺🇸 US FDA FCS record found for PubChem cid {cid}.")
            return None

    def __iter__(self):
        for rid in self.order:
            yield self._load_record(rid, order=rid)

    def __len__(self):
        return len(self.order)

    def __contains__(self, item):
        if isinstance(item, list):
            item = item[0]
        if isinstance(item, int):
            return item in self.order or (("bycid" in self.index) and (str(item) in self.index["bycid"]))
        if isinstance(item, str):
            return item in self.index.get("CAS", {}) or item in self.index.get("FCNNo", {})
        return False

    def __repr__(self):
        csv_filename = os.path.basename(self.csv_file)
        index_date = self.index.get("index_date", "unknown")
        print(f"🇺🇸US FDA FCS database ({len(self.order)} records)")
        print(f"Imported from CSV {csv_filename} and indexed on {index_date}")
        return str(self)

    def __str__(self):
        return f"<{self.__class__.__name__}: {len(self.order)} records (US FDA FCS)>"


# ----------------------------------------------------------------------
# Standalone test / debugging section
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # For debugging or standalone tests, one can initialize the database:
    db = USFDAfcn()
    print(db)
    # Example lookup:
    try:
        rec = db[1]
        print(rec)
    except Exception as e:
        print(f"Error retrieving record 1: {e}")
