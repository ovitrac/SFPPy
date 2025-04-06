"""
Module: GBappendixA

This module implements a robust database manager for the Chinese positive list under GB 9685-2016.
It is modeled after the EUFCMannex1 module but with the following key differences:

1. The CSV file (GB9685-2016.csv) merges 8 tables corresponding to different food contact materials.
   The table identifier is provided in the first column ("表格") and is mapped to descriptive names
   (e.g. "A1" becomes "plastics", "A2" becomes "coatings", etc.).
2. The primary key is based on the second column ("FCA编号"). From an entry like FCAXXXX the digits
   XXXX are extracted and stored under the key "FCA". Cached files are named FCAXXXX.json.
3. A substance may appear several times (once per material category) so that the positive list
   details (usage range, limits, restrictions, etc.) are stored in a sub-dictionary. The “authorized in”
   field is built as a list of the descriptive table names.
4. The parsing rules for the remaining columns are as follows:
   - Column 3 (“中文名称”): stored as "ChineseName".
   - Column 4 (“CAS号”): stored as a string or a list if multiple CAS numbers (separated by “;”).
   - Column 5 (“使用范围和最大使用量/%”): split into two fields within the positive list info:
       - "materials": a list obtained by splitting the text before “:” using commas.
       - "CP0max": the numeric value after “:” converted from a percentage (w/w) to mg/kg via
         multiplication by $10^4$. If no “:” is present, this field is set to None.
   - Column 6 (“SML/QM/(mg/kg)”): the full raw string is stored in "QMSMLraw". It is parsed to extract:
       - "SML": a single or list of numeric values (ignoring “ND” entries) when the associated text contains “:SML”.
       - "QM": a numeric value when the associated text contains “:QM”.
       - "DL": a detection limit if a pattern like “DL=0.01mg/kg” is found.
   - Column 7 (“SML(T)/(mg/kg)”): if the field is numeric (or a semicolon‑separated list of numbers) the
         numeric value(s) are stored in "SMLT"; otherwise, the raw text is stored in "SMLTraw".
   - Column 8 (“SML(T)”): stored as "SMLTcomment".
   - Column 9 (“分组编号”): stored as "comment1".
   - Column 10 (“其他要求”): stored as "comment2".
5. Finally, the order of the fields in each record is maintained as: (1) FCA, (2) cid, (3) CAS, (4) "authorized in",
   (5) ChineseName, then the positive list details (under "table"), and finally the traceability fields "engine", "csfile" and "date".

The extended record class, *gbrecord_ext*, automatically adds additional chemical information from PubChem.
The main manager class, *GBappendixA*, builds an index for fast lookup by CAS, FCA (the primary key), or PubChem cid.

@version: 1.41
@project: SFPPy - SafeFoodPackaging Portal in Python initiative
@author: INRAE\\olivier.vitrac@agroparistech.fr
@licence: MIT
@Date: 2024-01-10
@rev: 2025-04-01

"""

import os, csv, json, datetime, time, re, textwrap

__all__ = ['GBappendixA', 'custom_wrap', 'extract_number_before_keyword', 'extract_number_before_keyword_in_parentheses', 'gbrecord', 'gbrecord_ext', 'printWARN', 'split_col5_content', 'unwrap']


__project__ = "SFPPy"
__author__ = "Olivier Vitrac"
__copyright__ = "Copyright 2022"
__credits__ = ["Olivier Vitrac"]
__license__ = "MIT"
__maintainer__ = "Olivier Vitrac"
__email__ = "olivier.vitrac@agroparistech.fr"
__version__ = "1.41"


# Module-level variables to track last warning message and its timestamp
_LAST_WARN_ = None
_T_LAST_WARN_ = 0.0

# ----------------------------------------------------------------------
# Custom text wrapping function (similar to EU module)
def custom_wrap(text, width=60, indent=" " * 22):
    # Wrap the first line and indent subsequent lines.
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
# robust numeric extract for GB document
def unwrap(value):
    """
    Unwrap a list-like value:
    - Return None if the value is empty or None
    - Return the sole item if it is a singleton list or tuple
    - Return the full list/tuple otherwise
    """
    if not value:
        return None
    if isinstance(value, (list, tuple)) and len(value) == 1:
        return value[0]
    return value

def extract_number_before_keyword_in_parentheses(text, keyword):
    """
    Extract numbers that are immediately followed by parentheses containing a given keyword.

    Example:
        extract_number_before_keyword_in_parentheses(
            'this is 0.6 (SML) 0.5 (again SML)', 'SML'
        )
        returns: [0.6, 0.5]

    Behavior:
        - If one match is found: return the number as float
        - If no match is found: return None
        - If multiple matches: return list of floats
    """
    pattern = rf'([+-]?\d*\.?\d+)\s*\(([^)]*{re.escape(keyword)}[^)]*)\)'
    matches = re.findall(pattern, text, re.IGNORECASE)
    res = unwrap([float(num) for num, _ in matches])
    return extract_number_before_keyword(text, keyword) if res is None else res


def extract_number_before_keyword(text, keyword):
    """
    Extract numbers that immediately precede a given keyword in the text, allowing for optional punctuation
    or other separators between the number and the keyword (e.g. '3.4 (SML)', '3.4: SML', etc.)

    Example:
        extract_number_before_keyword('1.2 3.4 SML', 'SML')         → 3.4
        extract_number_before_keyword('1 SML 2.0, SML again', 'SML') → [1, 2.0]
        extract_number_before_keyword('4.5: SML and (5.5) SML', 'SML') → [4.5, 5.5]

    Returns:
        - A float or int if a single number is found
        - A list of floats/ints if multiple matches are found
        - None if no match
    """
    pattern = rf'([+-]?\d*\.?\d+)\s*[\W]*\s*{re.escape(keyword)}\b'
    matches = re.findall(pattern, text, re.IGNORECASE)
    return unwrap([float(m) if '.' in m else int(m) for m in matches])


def split_col5_content(text):
    """
    Pattern to match number followed by (...) or [...] that contains a target keyword
    """
    # True regulatory keywords
    keywords = r'(SML|DL|QM|SML\(T\))'
    # Match a value (ND or number) followed by (...) or [...] containing a keyword
    # We will match starting FROM that value
    pattern = re.compile(
        r'(ND|\d+(?:\.\d+)?)(\([^)]*' + keywords + r'[^)]*\)|\[[^\]]*' + keywords + r'[^\]]*\])'
    )
    match = pattern.search(text)
    if match:
        idx = match.start()  # this is the correct split point: just before the match
        main = text[:idx].strip()
        rem = text[idx:].strip()
        return main, rem
    else:
        return text.strip(), ""  # fallback


# ----------------------------------------------------------------------
# gbrecord: represents one substance record from GB 9685-2016.
# The record stores the following keys in order:
#   "FCA" (primary key extracted from FCA编号),
#   "cid" (PubChem identifier),
#   "CAS" (CAS number, string or list),
#   "authorized in" (list of material names in which the substance is allowed),
#   "ChineseName" (the substance name in Chinese),
#   "table" (a list of dictionaries with the positive list details),
#   and traceability keys "engine", "csfile", "date".
# ----------------------------------------------------------------------
class gbrecord(dict):
    def __init__(self, d, order=None, total=None):
        if not isinstance(d, dict):
            raise TypeError(f"dict must be a dict not a {type(d).__name__}")
        super().__init__(d)
        self._order = d.get("FCA", order)
        self._total = total

    def __str__(self):
        cid = self.get("cid", None)
        order_str = f"{self._order}" if self._order is not None else "?"
        total_str = f"{self._total}" if self._total is not None else "?"
        return f"<{self.__class__.__name__} with cid:{cid} - FCA {order_str} of {total_str} (GB 9685-2016)>"

    def __repr__(self):
        lines = []
        order_str = f"{self._order}" if self._order is not None else "?"
        total_str = f"{self._total}" if self._total is not None else "?"
        header = f" ---- [ GB 9685-2016 record: {order_str} of {total_str} ] ----"
        lines.append(header)

        # Define base keys to display first.
        base_keys = {"FCA", "cid", "CAS", "authorized in", "ChineseName", "engine", "csfile", "date"}
        fields_order = ["FCA", "cid", "CAS", "authorized in", "ChineseName"]
        for key in fields_order:
            if key not in self:
                continue
            val = self[key]
            if val is None or (isinstance(val, str) and not val.strip()):
                continue
            wrapped_val = custom_wrap(str(val), width=60, indent=" " * 22)
            lines.append(f"{key:>20}: {wrapped_val}")

        # Display branching positive list details:
        branch_keys = [k for k in self.keys() if k not in base_keys]
        if branch_keys:
            lines.append(f"\n{'--- Positive Lists':>20}:")
            for branch in sorted(branch_keys):
                value = self[branch]
                lines.append(f"  {branch}:")
                # Normalize value to a list
                entries = value if isinstance(value, list) else [value]
                for idx, entry in enumerate(entries, start=1):
                    lines.append(f"    Entry {idx}:")
                    if isinstance(entry, dict):
                        for k, v in entry.items():
                            if v is None or (isinstance(v, str) and not v.strip()):
                                continue
                            wrapped_v = custom_wrap(str(v), width=60, indent=" " * 25)
                            lines.append(f"      {k:>15}: {wrapped_v}")
                    else:
                        wrapped_entry = custom_wrap(str(entry), width=60, indent=" " * 25)
                        lines.append(f"      {wrapped_entry}")

        # Append traceability information.
        for key in ["engine", "csfile", "date"]:
            if key in self and self[key]:
                wrapped_val = custom_wrap(str(self[key]), width=60, indent=" " * 22)
                lines.append(f"{key:>20}: {wrapped_val}")
        return "\n".join(lines)


# ----------------------------------------------------------------------
# gbrecord_ext: extended record that adds chemical info from PubChem.
# ----------------------------------------------------------------------
class gbrecord_ext(gbrecord):
    def __init__(self, rec, db=None):
        if not isinstance(rec, gbrecord):
            raise TypeError(f"dict must be a gbrecord not a {type(rec).__name__}")
        super().__init__(rec, order=rec._order, total=rec._total)
        from patankar.loadpubchem import migrant
        if rec.get("CAS") and (isinstance(rec.get("CAS"), str) and rec.get("CAS").strip() or isinstance(rec.get("CAS"), list)):
            try:
                m = migrant(rec.get("CAS"), annex1=False)
                M = m.M
            except Exception:
                print(f"{rec.get('ChineseName')} could not be extended from its CAS {rec.get('CAS')}: not found")
                M = None
        else:
            M = None
        self.cid = self.get("cid")
        self.M = M
        # Additional extended attributes (e.g. group information) can be added here.
        # [XXX] Additional placeholders remain as needed.

# ----------------------------------------------------------------------
# GBappendixA: the main class to manage the GB 9685-2016 CSV file and cache.
#
# The parsing of each CSV row follows these rules:
#
# 1) Column 1: 表格
#    - The code (e.g. "A1") is mapped to a descriptive name (e.g. "plastics") and stored in the record under "authorized in".
#
# 2) Column 2: FCA编号
#    - Expected in the format "FCAXXXX"; the digits XXXX are extracted and stored under "FCA".
#
# 3) Column 3: 中文名称
#    - Stored as a string under "ChineseName".
#
# 4) Column 4: CAS号
#    - Stored as a string or as a list (if multiple numbers are provided separated by ";").
#
# 5) Column 5: 使用范围和最大使用量/%
#    - Parsed to extract:
#         * "materials": a list (splitting by commas before the colon).
#         * "CP0max": the numeric value after the colon converted to mg/kg by multiplying by $10^4$.
#       If no colon is present, "CP0max" is set to None.
#
# 6) Column 6: SML/QM/(mg/kg)
#    - The full text is stored in "QMSMLraw".
#    - It is parsed by splitting on ";" to detect multiple entries.
#      For each entry:
#         * If the entry contains ":SML" and the value is not "ND", the number is added to the SML list.
#         * If the entry contains ":QM", the numeric value is recorded as QM.
#         * If a pattern like "DL=0.01mg/kg" is found, DL is extracted.
#
# 7) Column 7: SML(T)/(mg/kg)
#    - If the field can be converted to a float (or is a semicolon‑separated list of floats), the value(s) are stored in "SMLT";
#      otherwise, the raw text is stored in "SMLTraw".
#
# 8) Column 8: SML(T)
#    - Stored as "SMLTcomment".
#
# 9) Column 9: 分组编号
#    - Stored as "comment1".
#
# 10) Column 10: 其他要求
#     - Stored as "comment2".
#
# If a substance appears in several rows (i.e. in different tables), the record is merged:
#    - The "authorized in" field becomes the list of all descriptive table names.
#    - The positive list details (columns 5–10) are stored as a list under "table".
#
# Traceability fields "engine", "csfile" and "date" are appended to each record.
#
# A global index is built mapping keys for "CAS", "FCA", "bycid" (PubChem cid) and "ChineseName" to the corresponding record.
# ----------------------------------------------------------------------
class GBappendixA:
    def __init__(self, cache_dir="cache.GBappendixA", index_file="gb_index.json", pubchem=True):
        self.base_dir = os.path.dirname(__file__)
        self.csv_file = os.path.join(self.base_dir, "GB9685-2016.csv")
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
        self._pubchem = pubchem # we enforce pubchem, the database is initialized indeed
        GBappendixA.isinitialized = True

    @classmethod
    def isindexinitialized(cls, cache_dir="cache.GBappendixA", index_file="gb_index.json"):
        return os.path.exists(os.path.join(os.path.dirname(__file__), cache_dir, index_file))

    def refresh_index(self):
        """
        Rebuild the global index by reading the CSV file and regenerating each record as FCAXXXX.json.
        The index includes mappings for "CAS", "FCA", "bycid", and "ChineseName".
        """
        # Load missing CAS numbers (only for substances with a CAS)
        missing_file = os.path.join(self.cache_dir, "missing.pubchem.gb.json")
        if os.path.exists(missing_file):
            with open(missing_file, "r", encoding="utf-8") as mf:
                missing_pubchem = json.load(mf)
        else:
            missing_pubchem = {}

        new_index = {}
        new_index["index_date"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        new_index["csv_file"] = os.path.basename(self.csv_file)
        new_index["order"] = []
        new_index["CAS"] = {}
        new_index["bycid"] = {}
        new_index["FCA"] = {}
        new_index["ChineseName"] = {}

        # Temporary dictionary to merge records by FCA number.
        records_dict = {}
        # Mapping table codes to descriptive names.
        table_mapping = {
            "A1": "plastics",
            "A2": "coatings",
            "A3": "rubber",
            "A4": "printing inks",
            "A5": "adhesives",
            "A6": "paper and board",
            "A7": "silicon rubber",
            "A7bis": "textile"
        }

        from patankar.loadpubchem import migrant  # Import the PubChem lookup function

        with open(self.csv_file, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter=",", quotechar='"')
            header = next(reader, None)
            # Assume header row exists (starting with "表格")
            for row in reader:
                if not row or len(row) < 10:
                    continue

                # Column 1: 表格
                table_code = row[0].strip()
                table_desc = table_mapping.get(table_code, table_code)

                # Column 2: FCA编号 (e.g. "FCA0001")
                fca_field = row[1].strip()
                m_fca = re.search(r'FCA(\d+)', fca_field)
                if m_fca:
                    fca_num = m_fca.group(1)
                else:
                    continue  # Skip row if FCA format is not recognized.

                # Column 3: 中文名称
                chinese_name = row[2].strip()

                # Column 4: CAS号
                cas_field = row[3].strip()
                if ";" in cas_field:
                    cas_value = [x.strip() for x in cas_field.split(";") if x.strip()]
                else:
                    cas_value = cas_field

                # Column 5: 使用范围和最大使用量/%
                range_field,remcol5 = split_col5_content(row[4].strip())
                materials = []
                CP0max = None
                if ":" in range_field:
                    parts = range_field.split(":")
                    materials = [x.strip() for x in parts[0].split(",") if x.strip()]
                    try:
                        CP0max = float(parts[1].strip()) * 1e4
                    except:
                        CP0max = None
                else:
                    materials = [x.strip() for x in range_field.split(",") if x.strip()]

                # Column 6: SML/QM/(mg/kg)
                smlqm_field = remcol5 + row[5].strip()
                QMSMLraw = smlqm_field
                sml_values = []
                qm_value = None
                dl_value = None
                entries = [entry.strip() for entry in smlqm_field.split(";") if entry.strip()]
                for entry in entries:
                    entry = entry.rstrip("或")
                    m_entry = re.match(r'(?P<val>\d*\.?\d+|ND)\s*\((?P<info>.+?)\)', entry)
                    if m_entry:
                        val_str = m_entry.group("val")
                        info = m_entry.group("info")
                        if val_str.upper() != "ND":
                            try:
                                num_val = float(val_str)
                            except:
                                num_val = None
                        else:
                            num_val = None
                        if ":SML" in info:
                            if num_val is not None:
                                sml_values.append(num_val)
                        if ":QM" in info:
                            if num_val is not None:
                                qm_value = num_val
                        dl_match = re.search(r'DL=([\d\.]+)mg/kg', info)
                        if dl_match:
                            try:
                                dl_value = float(dl_match.group(1))
                            except:
                                dl_value = None
                SML = unwrap(sml_values)
                QM = unwrap(qm_value)
                DL = unwrap(dl_value)
                # second pass for column 6
                if SML is None and "SML" in smlqm_field:
                    SML = extract_number_before_keyword_in_parentheses(smlqm_field,"SML")
                if QM is None and "QM" in smlqm_field:
                    QM = extract_number_before_keyword_in_parentheses(smlqm_field, "QM")
                if DL is None and "DL" in smlqm_field:
                    DL = extract_number_before_keyword_in_parentheses(smlqm_field, "DL")

                # Column 7: SML(T)/(mg/kg)
                smlt_field = row[6].strip()
                SMLT = None
                SMLTraw = ""
                if smlt_field:
                    parts = [x.strip() for x in smlt_field.split(";") if x.strip()]
                    parsed_parts = []
                    for part in parts:
                        try:
                            parsed_parts.append(float(part))
                        except:
                            parsed_parts.append(part)
                    if all(isinstance(x, float) for x in parsed_parts):
                        SMLT = parsed_parts[0] if len(parsed_parts) == 1 else parsed_parts
                    else:
                        SMLT = extract_number_before_keyword_in_parentheses(smlt_field, "SML")
                        SMLTraw = smlt_field if SMLT is None else SMLTraw

                # consistency rule
                if SMLT is not None and SML is None:
                    SML = SMLT

                # Column 8: SML(T)
                SMLTcomment = row[7].strip()
                # Column 9: 分组编号
                comment1 = row[8].strip()
                # Column 10: 其他要求
                comment2 = row[9].strip()

                # Assemble the positive list info for this row.
                pos_info = {
                    "materials": materials,
                    "CP0max": CP0max,
                    "QMSMLraw": QMSMLraw,
                    "SML": SML,
                    "QM": QM,
                    "DL": DL,
                    "SMLT": SMLT,
                    "SMLTraw": SMLTraw,
                    "SMLTcomment": SMLTcomment,
                    "comment1": comment1,
                    "comment2": comment2,
                    "table_id": table_desc
                }

                # --- Begin PubChem lookup for cid ---
                if cas_value:
                    if isinstance(cas_value, list):
                        cas_lookup = cas_value[0]
                    else:
                        cas_lookup = cas_value
                    if cas_lookup and cas_lookup.strip():
                        if cas_lookup in missing_pubchem:
                            cid_val = missing_pubchem[cas_lookup]
                        else:
                            try:
                                cid_val = migrant(cas_lookup, annex1=False).cid
                            except ValueError:
                                printWARN(f"🇨🇳 Warning: substance {chinese_name} (CAS {cas_lookup}) not found in PubChem.")
                                cid_val = None
                                missing_pubchem[cas_lookup] = None
                    else:
                        cid_val = None
                else:
                    cid_val = None
                # --- End PubChem lookup ---

                # Prepare a temporary record for this row using columns 1-4 and traceability.
                temp_rec = {
                    "FCA": fca_num,
                    "cid": cid_val,
                    "CAS": cas_value,
                    "authorized in": [table_desc],
                    "ChineseName": chinese_name,
                    "engine": "SFPPy: GBappendixA module",
                    "csfile": os.path.basename(self.csv_file),
                    "date": new_index["index_date"]
                }
                # The positive list details (columns 5-10) are stored in pos_info.
                # Merge records by FCA number.
                if fca_num in records_dict:
                    record = records_dict[fca_num]
                    # Update the "authorized in" list if the category is new.
                    if table_desc not in record.get("authorized in", []):
                        record["authorized in"].append(table_desc)
                    # Branch under the key corresponding to table_desc.
                    if table_desc in record:
                        # If an entry already exists for this category, append the new pos_info.
                        existing = record[table_desc]
                        if isinstance(existing, list):
                            existing.append(pos_info)
                        else:
                            record[table_desc] = [existing, pos_info]
                    else:
                        record[table_desc] = pos_info
                else:
                    # For a new FCA, include the branch with key table_desc.
                    temp_rec[table_desc] = pos_info
                    records_dict[fca_num] = temp_rec

                # Update index for CAS and ChineseName.
                if cas_value:
                    if isinstance(cas_value, list):
                        for cas in cas_value:
                            new_index["CAS"].setdefault(cas, []).append(fca_num)
                    else:
                        new_index["CAS"].setdefault(cas_value, []).append(fca_num)
                new_index["FCA"].setdefault(fca_num, []).append(fca_num)
                new_index["ChineseName"].setdefault(chinese_name, []).append(fca_num)

        # Write individual record files and build the order list.
        order_list = sorted(records_dict.keys(), key=lambda x: int(x))
        for fca in order_list:
            record = records_dict[fca]
            record_filename = f"FCA{int(fca):04d}.json"
            json_filename = os.path.join(self.cache_dir, record_filename)
            with open(json_filename, "w", encoding="utf-8") as jf:
                json.dump(record, jf, ensure_ascii=False, indent=2)
            new_index["order"].append(fca)
            if record.get("cid") is not None:
                new_index["bycid"][str(record["cid"])] = fca

        with open(self.index_file, "w", encoding="utf-8") as f:
            json.dump(new_index, f, ensure_ascii=False, indent=2)
        with open(missing_file, "w", encoding="utf-8") as mf:
            json.dump(missing_pubchem, mf, ensure_ascii=False, indent=2)
        self.index = new_index
        self.order = new_index.get("order", [])
        self._records_cache = {}


    def _load_record(self, fca, order=None, db=False):
        """
        Load a record (as a gbrecord) from its cached JSON file.
        If PubChem extension is enabled, the record is returned as a gbrecord_ext.
        """
        if fca in self._records_cache:
            if self._pubchem:
                return gbrecord_ext(self._records_cache[fca], self) if db else gbrecord_ext(self._records_cache[fca])
            else:
                return self._records_cache[fca]
        json_filename = os.path.join(self.cache_dir, f"FCA{int(fca):04d}.json")
        if not os.path.exists(json_filename):
            print(f"⚠️ Warning: Record file for 🇨🇳 FCA {fca} not found.")
            return None
        with open(json_filename, "r", encoding="utf-8") as jf:
            rec = json.load(jf)
        record_obj = gbrecord(rec, order=rec.get("FCA"), total=len(self.order))
        self._records_cache[fca] = record_obj
        if self._pubchem:
            return gbrecord_ext(record_obj, self) if db else gbrecord_ext(record_obj)
        else:
            return record_obj

    def __getitem__(self, key):
        """
        Supports lookup by:
         - Slice: returns a list of records whose FCA numbers fall within the slice.
         - Integer or string: interpreted as an FCA number.
         - String: if not matching an FCA number then interpreted as a CAS number.
         - List/tuple: returns a list of corresponding records.
        """
        if isinstance(key, slice):
            start = key.start if key.start is not None else min(self.order, key=lambda x: int(x))
            stop = key.stop if key.stop is not None else max(self.order, key=lambda x: int(x)) + 1
            rec_keys = [k for k in self.order if int(start) <= int(k) < int(stop)]
            if not rec_keys:
                raise KeyError(f"No records found in range {start} to {stop - 1}. Valid FCA numbers range from {min(self.order)} to {max(self.order)}.")
            return [self._load_record(k, order=k) for k in rec_keys]
        elif isinstance(key, (int, str)):
            key_str = str(key)
            if key_str in self.order:
                return self._load_record(key_str, order=key_str)
            elif key in self.index.get("CAS", {}):
                rec_keys = self.index["CAS"][key]
                if len(rec_keys) == 1:
                    return self._load_record(rec_keys[0], order=rec_keys[0])
                else:
                    return [self._load_record(k, order=k) for k in rec_keys]
            else:
                raise KeyError(f"Key '{key}' not found. Valid keys include FCA numbers and CAS numbers.")
        elif isinstance(key, (list, tuple)):
            return [self.__getitem__(k) for k in key]
        else:
            raise KeyError(f"Unsupported key type: {type(key)}")

    def __call__(self, *args):
        """
        Callable access. For example:
         - GBappendixA(cid) returns the record for a given PubChem cid.
         - GBappendixA(fca) returns the record for a given FCA number.
         - GBappendixA("CAS") returns the record(s) for a given CAS number.
         - Multiple arguments return a list.
        """
        if len(args) == 1 and isinstance(args[0], (list, tuple)):
            args = args[0]
        results = []
        for arg in args:
            if isinstance(arg, int):
                arg_str = str(arg)
                if arg_str in self.order:
                    results.append(self._load_record(arg_str))
                elif arg_str in self.index.get("bycid", {}):
                    fca = self.index["bycid"][arg_str]
                    results.append(self._load_record(fca))
                else:
                    print(f"🇨🇳 Warning: Record for identifier {arg} not found.")
                    results.append(None)
            elif isinstance(arg, str):
                try:
                    result_item = self.__getitem__(arg)
                    if isinstance(result_item, list):
                        results.extend(result_item)
                    else:
                        results.append(result_item)
                except KeyError as e:
                    print(e)
                    results.append(None)
            else:
                raise KeyError(f"Unsupported key type in call: {type(arg)}")
        return results[0] if len(results) == 1 else results

    def byCAS(self, cas):
        if isinstance(cas, list):
            cas = cas[0]
        rec_keys = self.index.get("CAS", {}).get(cas, [])
        if len(rec_keys) == 1:
            return self._load_record(rec_keys[0], order=rec_keys[0], db=True)
        else:
            return [self._load_record(k, order=k, db=True) for k in rec_keys]

    def byFCA(self, fca):
        fca_str = str(fca)
        if fca_str in self.order:
            return self._load_record(fca_str, order=fca_str)
        else:
            raise KeyError(f"🇨🇳 FCA number {fca} not found. Valid FCA numbers range from {min(self.order)} to {max(self.order)}.")

    def bycid(self, cid, verbose=True):
        cid_str = str(cid)
        if cid_str in self.index.get("bycid", {}):
            fca = self.index["bycid"][cid_str]
            return self._load_record(fca, order=fca, db=True)
        else:
            if verbose:
                print(f"⚠️ Warning: No 🇨🇳 GB 9685-2016 record found for PubChem cid {cid}.")
            return None

    def __iter__(self):
        for fca in self.order:
            yield self._load_record(fca, order=fca)

    def __len__(self):
        return len(self.order)

    def __contains__(self, item):
        if isinstance(item, (list, tuple)):
            item = item[0]
        if isinstance(item, int):
            return str(item) in self.order or str(item) in self.index.get("bycid", {})
        if isinstance(item, str):
            return item in self.index.get("CAS", {})
        return False

    def __repr__(self):
        csv_filename = os.path.basename(self.csv_file)
        index_date = self.index.get("index_date", "unknown")
        print(f"GB 9685-2016 positive list ({len(self.order)} records)")
        print(f"Imported from CSV {csv_filename} and indexed on {index_date}")
        return str(self)

    def __str__(self):
        return f"<{self.__class__.__name__}: {len(self.order)} records (GB 9685-2016)>"

# -------------------------------------------------------------------
# Example usage (for debugging / standalone tests)
# -------------------------------------------------------------------
if __name__ == "__main__":
    dbappendixA = GBappendixA(pubchem=True)
    # Lookup by FCA number:
    rec = dbappendixA.byFCA("0001")
    print(repr(rec))
    # Lookup by CAS:
    try:
        rec_by_cas = dbappendixA.byCAS("25013-16-5")
        print(rec_by_cas)
    except KeyError as e:
        print(e)
    # Lookup by PubChem cid:
    rec_by_cid = dbappendixA.bycid(6581)
    print(rec_by_cid)
    # Callable access example:
    rec_call = dbappendixA("25013-16-5")
    print(rec_call)
    # Iterate over records:
    #for record in dbappendixA:
    #    print(record.get("FCA"), record.get("ChineseName"))
