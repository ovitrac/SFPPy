"""
Substances Routes - Chemical Substance Management

API endpoints for substance search, PubChem integration, and C0/SML management.
"""

import sys
from pathlib import Path
from typing import List, Optional, Dict, Any

from fastapi import APIRouter, HTTPException, Query, Response
from fastapi.responses import JSONResponse, FileResponse
from pydantic import BaseModel, Field
import os

# Add parent paths
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

router = APIRouter()


# ========== HELPERS ==========

def _extract_scalar(value):
    """Extract scalar from array/list/numpy array, or return the value itself."""
    if value is None:
        return None
    # Handle numpy arrays
    if hasattr(value, '__iter__') and not isinstance(value, str):
        try:
            return float(value[0]) if len(value) > 0 else None
        except (IndexError, TypeError):
            return None
    return float(value) if value is not None else None


# ========== DATA MODELS ==========

class SubstanceInput(BaseModel):
    name: str = Field(..., description="Substance name or identifier")
    cas: Optional[str] = Field(None, description="CAS registry number")
    cid: Optional[int] = Field(None, description="PubChem CID")
    mw: Optional[float] = Field(None, description="Molecular weight (g/mol)")
    logP: Optional[float] = Field(None, description="Log partition coefficient")
    smiles: Optional[str] = Field(None, description="SMILES structure")


class SubstanceInLayer(BaseModel):
    substance_id: str = Field(..., description="Substance identifier")
    layer_index: int = Field(..., ge=1, le=10, description="Layer containing substance")
    C0: float = Field(..., gt=0, description="Initial concentration (mg/kg)")
    SML: Optional[float] = Field(None, description="Specific migration limit (mg/kg food)")
    use_SML_check: bool = Field(default=True, description="Enable SML compliance check")


class SubstanceSearchResult(BaseModel):
    cid: int
    name: str
    iupac_name: Optional[str]
    mw: float
    formula: str
    smiles: Optional[str]
    logP: Optional[float]
    thumbnail_url: Optional[str]


# ========== PUBCHEM INTEGRATION ==========

def get_pubchem_loader():
    """Get PubChem loader via migrant class."""
    try:
        from patankar.loadpubchem import migrant
        return migrant
    except ImportError:
        return None


def get_migrant_class():
    """Get migrant class for basic substance info."""
    try:
        from patankar.loadpubchem import migrant
        return migrant
    except ImportError:
        return None


def get_migrant_toxtree_class():
    """Get migrantToxtree class for full toxicological data."""
    try:
        from patankar.loadpubchem import migrantToxtree
        return migrantToxtree
    except ImportError:
        return None


def load_substance_with_toxtree(query: str) -> Optional[Dict[str, Any]]:
    """
    Load substance with full toxicological and regulatory data from migrantToxtree.

    Returns comprehensive data including:
    - PubChem properties (name, CAS, MW, formula, SMILES, logP, InChiKey)
    - EU 10/2011 regulatory status (SML, FCM No, REF, authorized)
    - US FCN list data (FCM No, Notifier, Manufacturer)
    - CN GB9685-2016 data (FCA No, authorized_in, polymers, SML)
    - Cramer classification (class, TTC, CF_TTC)
    - Toxicological structural alerts

    Handles both name-based queries and CID lookups (numeric strings).
    """
    toxtree_cls = get_migrant_toxtree_class()
    migrant_cls = get_migrant_class()

    if toxtree_cls is None and migrant_cls is None:
        return None

    original_query = query
    mt = None

    try:
        # Try migrantToxtree first for full data
        if toxtree_cls is not None:
            try:
                mt = toxtree_cls(query, verbose=False)
            except Exception:
                pass

        # Fallback to basic migrant if toxtree failed
        if mt is None and migrant_cls is not None:
            try:
                m = migrant_cls(query, verbose=False)
                if m is not None and hasattr(m, 'cid') and m.cid:
                    # Try to promote to toxtree
                    if hasattr(m, 'ispromovable') and m.ispromovable() and hasattr(m, 'promote'):
                        mt = m.promote()
                    else:
                        # Use basic migrant data
                        return _build_basic_result(m, original_query)
            except Exception:
                pass

        if mt is None or not hasattr(mt, 'cid') or mt.cid is None:
            return None

        # Extract the best name - prefer user's query if it matches
        name = _get_best_name(mt, original_query)

        # Get a few good synonyms for display
        synonyms = _get_synonyms(mt, name, max_synonyms=3)

        # Extract CAS (can be string or list)
        cas = None
        if hasattr(mt, 'CAS') and mt.CAS:
            cas = mt.CAS if isinstance(mt.CAS, list) else [mt.CAS]
        elif hasattr(mt, 'cas') and mt.cas:
            cas = mt.cas if isinstance(mt.cas, list) else [mt.cas]

        # Build comprehensive result
        result = {
            "cid": mt.cid,
            "name": name,
            "synonyms": synonyms,
            "iupac_name": mt.iupac_name if hasattr(mt, 'iupac_name') else None,
            "cas": cas,
            "mw": mt.M if hasattr(mt, 'M') else None,
            "formula": mt.formula if hasattr(mt, 'formula') else None,
            "smiles": mt.smiles if hasattr(mt, 'smiles') else None,
            "inchikey": mt.InChiKey if hasattr(mt, 'InChiKey') else None,
            "logP": _extract_scalar(mt.logP) if hasattr(mt, 'logP') else None,
            "thumbnail_url": f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={mt.cid}&t=l",
            "structure_2d_url": f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={mt.cid}&t=s",
        }

        # EU 10/2011 regulatory data
        # SML is a top-level attribute, hasSML indicates if it's set
        eu_data = {"authorized": None, "SML": None, "FCM": None, "REF": None, "EC": None}
        has_sml = getattr(mt, 'hasSML', False)
        sml_value = getattr(mt, 'SML', None)
        if has_sml and sml_value is not None:
            eu_data["authorized"] = True
            eu_data["SML"] = float(sml_value)
        # Try specific EU* attributes too
        if hasattr(mt, 'EUauthorized') and mt.EUauthorized:
            eu_data["authorized"] = True
        if hasattr(mt, 'EUSML') and mt.EUSML:
            eu_data["SML"] = float(mt.EUSML)
        if hasattr(mt, 'EUFCMNo') and mt.EUFCMNo:
            eu_data["FCM"] = mt.EUFCMNo
        if hasattr(mt, 'EUREFNo') and mt.EUREFNo:
            eu_data["REF"] = mt.EUREFNo
        if hasattr(mt, 'EUECNo') and mt.EUECNo:
            eu_data["EC"] = mt.EUECNo
        result["regulatory_EU"] = eu_data

        # US FCN list data
        # FCNNo is a top-level attribute (list of FCM numbers)
        us_data = {"authorized": None, "FCM_No": None, "Notifier": None, "Manufacturer": None, "N_Date": None, "Mixture": None}
        has_fcn = getattr(mt, 'hasFCN', False)
        fcn_no = getattr(mt, 'FCNNo', None)
        if has_fcn and fcn_no:
            us_data["authorized"] = True
            us_data["FCM_No"] = fcn_no
        # Try specific US* attributes too
        if hasattr(mt, 'USauthorized') and mt.USauthorized:
            us_data["authorized"] = True
        if hasattr(mt, 'USFCMNo') and mt.USFCMNo:
            us_data["FCM_No"] = mt.USFCMNo
        if hasattr(mt, 'USNotifier') and mt.USNotifier:
            us_data["Notifier"] = mt.USNotifier
        if hasattr(mt, 'USManufacturer') and mt.USManufacturer:
            us_data["Manufacturer"] = mt.USManufacturer
        if hasattr(mt, 'USNDate') and mt.USNDate:
            us_data["N_Date"] = mt.USNDate
        if hasattr(mt, 'USMixture') and mt.USMixture:
            us_data["Mixture"] = mt.USMixture
        result["regulatory_US"] = us_data

        # CN GB9685-2016 data
        # FCANo is a top-level attribute
        cn_data = {"authorized": None, "FCA_No": None, "authorized_in": None, "polymers": None, "SML": None}
        has_fca = getattr(mt, 'hasFCA', False)
        fca_no = getattr(mt, 'FCANo', None)
        if has_fca and fca_no:
            cn_data["authorized"] = True
            cn_data["FCA_No"] = fca_no
        # Try specific CN* attributes too
        if hasattr(mt, 'CNauthorized') and mt.CNauthorized:
            cn_data["authorized"] = True
        if hasattr(mt, 'CNFCANo') and mt.CNFCANo:
            cn_data["FCA_No"] = mt.CNFCANo
        if hasattr(mt, 'CNauthorizedIn') and mt.CNauthorizedIn:
            cn_data["authorized_in"] = mt.CNauthorizedIn
        if hasattr(mt, 'CNpolymers') and mt.CNpolymers:
            cn_data["polymers"] = mt.CNpolymers
        cn_sml = getattr(mt, 'CNSML', None)
        if cn_sml:
            cn_data["SML"] = float(cn_sml)
        result["regulatory_CN"] = cn_data

        # Combined regulatory object for backwards compatibility
        result["regulatory"] = {
            "EU": {
                "authorized": eu_data["authorized"],
                "SML": eu_data["SML"],
                "FCM": eu_data["FCM"],
                "REF": eu_data["REF"],
            },
            "US": {
                "authorized": us_data["authorized"],
                "FCM_No": us_data["FCM_No"],
            },
            "CN": {
                "authorized": cn_data["authorized"],
                "SML": cn_data["SML"],
                "FCA_No": cn_data["FCA_No"],
            },
        }

        # ToxTree data (Cramer classification, TTC)
        # These are direct attributes on migrantToxtree, not in the ToxTree dict
        toxtree_data = None

        cramer_class = getattr(mt, 'CramerClass', None)
        cramer_value = getattr(mt, 'CramerValue', None)
        ttc = getattr(mt, 'TTC', None)
        cf_ttc = getattr(mt, 'CFTTC', None)

        # Also try from ToxTree dict if direct attributes are None
        if hasattr(mt, 'ToxTree') and mt.ToxTree and isinstance(mt.ToxTree, dict):
            tt = mt.ToxTree
            if cramer_class is None:
                cramer_class = tt.get("CramerRules") or tt.get("CramerClass")
            if cramer_value is None:
                cramer_value = tt.get("CramerValue", 0)

        # Convert cramer value to integer if needed
        if isinstance(cramer_value, str):
            cramer_value = {"I": 1, "II": 2, "III": 3}.get(cramer_value, 0)

        # If we have any ToxTree data, build the result
        if cramer_class is not None or ttc is not None:
            # Compute TTC if not available (fallback values)
            # Class I (low concern) = 30 µg/kg bw/day
            # Class II (intermediate) = 9 µg/kg bw/day
            # Class III (high concern) = 1.5 µg/kg bw/day
            if ttc is None and cramer_value is not None:
                TTC_VALUES = {0: 0.0025, 1: 30.0, 2: 9.0, 3: 1.5}
                ttc = TTC_VALUES.get(cramer_value, 0.0025)

            # Compute CF_TTC if not available
            if cf_ttc is None and ttc is not None:
                cf_ttc = ttc * 60 * 1e-3  # mg/kg food (60 kg bw, 1 kg food/day)

            # Get alerts info from skin sensitization plugin
            has_alerts = getattr(mt, 'has_alerts', False)
            nalerts = getattr(mt, 'nalerts', 0)
            showalerts = getattr(mt, 'showalerts', {})

            # Parse alerts and check for genotoxic patterns
            alerts_list = []
            has_genotoxic = False
            genotox_patterns = ['genotox', 'mutagen', 'carcinogen', 'dna', 'michael']

            for alert_key, alert_text in showalerts.items():
                alert_text_lower = str(alert_text).lower()
                is_genotoxic = any(pattern in alert_text_lower for pattern in genotox_patterns)
                if is_genotoxic:
                    has_genotoxic = True
                alerts_list.append({
                    "id": alert_key,
                    "text": alert_text,
                    "is_genotoxic": is_genotoxic,
                    "source": "skin"
                })

            # Also check DNA binding plugin for genotoxicity if available
            try:
                if hasattr(mt, '_run_toxtree'):
                    dna_result = mt._run_toxtree('dnabinding')
                    if dna_result:
                        # Check for DNA binding alerts
                        dna_alert_keys = [
                            ('AlertForMichaelAcceptorIdentified', 'Michael Acceptor (DNA reactive)'),
                            ('AlertForAcylTransferAgentIdentified', 'Acyl Transfer Agent (DNA reactive)'),
                            ('AlertForSn1Identified', 'SN1 mechanism (DNA reactive)'),
                            ('AlertForSn2Identified', 'SN2 mechanism (DNA reactive)'),
                            ('AlertForSchiffBaseFormationIdentified', 'Schiff Base Formation (DNA reactive)'),
                        ]
                        for key, description in dna_alert_keys:
                            if dna_result.get(key, '').upper() == 'YES':
                                has_genotoxic = True
                                nalerts += 1
                                alerts_list.append({
                                    "id": f"DNA_{key}",
                                    "text": description,
                                    "is_genotoxic": True,
                                    "source": "dnabinding"
                                })
            except Exception:
                pass  # DNA binding plugin may not be available

            # Update has_alerts if DNA binding found alerts
            has_alerts = has_alerts or len(alerts_list) > 0

            # If genotoxic, apply stricter TTC threshold
            genotox_ttc = 0.0025  # µg/kg bw/day for genotoxic compounds
            genotox_cf_ttc = genotox_ttc * 60 * 1e-3  # mg/kg food

            toxtree_data = {
                "cramer_class": cramer_class,
                "cramer_value": cramer_value,
                "ttc": ttc,
                "ttc_unit": "µg/kg bw/day",
                "cf_ttc": cf_ttc,
                "cf_ttc_unit": "mg/kg food",
                "has_alerts": has_alerts,
                "nalerts": nalerts,
                "alerts": alerts_list,
                "has_genotoxic": has_genotoxic,
                "genotox_ttc": genotox_ttc if has_genotoxic else None,
                "genotox_cf_ttc": genotox_cf_ttc if has_genotoxic else None,
            }
        result["toxtree"] = toxtree_data

        # Also set top-level SML from EU data if available
        if eu_data["SML"]:
            result["SML"] = eu_data["SML"]
        elif cn_data["SML"]:
            result["SML"] = cn_data["SML"]
        else:
            result["SML"] = None

        return result

    except Exception as e:
        import traceback
        traceback.print_exc()
        return None


def _name_score(n: str) -> float:
    """Score a name - lower is better."""
    score = 0
    # Strong penalties for code-like names
    if n.startswith(('DTXSID', 'SCHEMBL', 'CHEBI', 'NSC', 'EINECS', 'EC ', 'UNII-', 'AKOS', 'RefChem', 'HMS', 'CCG-', 'NCGC', 'BRD-', 'STK', 'SR-', 'AB', 'DB-', 'FS-', 'HY-', 'CS-', 'LS-', 'EN300', 'FO', 'FA', 'DA-', 'NS0', 'STL', 'WLN:')):
        score += 200
    # Penalty for CAS-like patterns (digits with dashes)
    if any(c.isdigit() for c in n[:3]):
        score += 100
    # Penalty for names with brackets/annotations
    if '(' in n or '[' in n:
        score += 50
    # Penalty for ALL-CAPS codes
    if n.isupper() and len(n) > 5:
        score += 80
    # Penalty for very long names
    if len(n) > 30:
        score += len(n) - 30
    # Penalty for technical/commercial terms
    lower = n.lower()
    if any(w in lower for w in ['standard', 'impurity', 'reference', 'certified', 'analytical', 'reagent', 'solution', 'microg/ml', 'puriss', 'grade']):
        score += 40
    # Slight penalty for very short names (might be abbreviations)
    if len(n) < 3:
        score += 20
    # Tiny penalty for length to prefer shorter names when equal
    score += len(n) * 0.01
    return score


def _get_best_name(obj, original_query: str) -> str:
    """
    Extract the best name from a migrant object.

    Priority:
    1. If original query matches a synonym exactly (case-insensitive), use it
    2. If original query is a partial match of a synonym, use that synonym
    3. Otherwise, find the best recognizable chemical name
    """
    name = obj.name if hasattr(obj, 'name') else original_query
    if isinstance(name, list):
        query_lower = original_query.lower().strip()

        # First check for exact match (case-insensitive)
        for n in name:
            if n.lower() == query_lower:
                return n

        # Check for partial matches - if query is contained in a name or vice versa
        for n in name:
            n_lower = n.lower()
            # Query contained in name (e.g., "benzene" in "Benzene, AR, >=99.5%")
            if query_lower in n_lower:
                # But only if it's a recognizable word boundary
                if n_lower.startswith(query_lower) or f' {query_lower}' in n_lower:
                    return n.split(',')[0].strip()  # Return the base name part

        # Check if any name contains the query as a main term
        for n in name:
            n_lower = n.lower()
            # E.g., user searched "benzene", there's a "Benzene" somewhere
            if query_lower in n_lower and len(query_lower) >= 4:
                # Get just the matched part cleaned up
                return original_query.capitalize() if len(original_query) <= 15 else n.split(',')[0].strip()

        # Find the best name: prefer recognizable chemical names
        sorted_names = sorted(name, key=_name_score)
        return sorted_names[0] if sorted_names else original_query
    return name


def _get_synonyms(obj, primary_name: str, max_synonyms: int = 3) -> List[str]:
    """
    Get a few good synonyms for display, excluding the primary name.

    Returns recognizable chemical names, avoiding abbreviations and codes.
    """
    names = obj.name if hasattr(obj, 'name') else []
    if not isinstance(names, list):
        return []

    # Filter out the primary name
    primary_lower = primary_name.lower()
    candidates = [n for n in names if n.lower() != primary_lower]

    # Custom scoring for synonyms - prefer meaningful names, not abbreviations
    def synonym_score(n: str) -> float:
        score = 0
        # Penalize very short names (abbreviations like BNZ, Ph-H)
        if len(n) <= 5:
            score += 100
        # Penalize ALL-CAPS abbreviations
        if n.isupper():
            score += 80
        # Penalize codes/database IDs
        if n.startswith(('DTXSID', 'SCHEMBL', 'CHEBI', 'NSC', 'EINECS', 'UNII', 'AKOS')):
            score += 200
        # Penalize CAS-like patterns
        if any(c.isdigit() for c in n[:3]):
            score += 150
        # Prefer medium-length names (8-20 chars)
        if 8 <= len(n) <= 20:
            score -= 20
        # Penalize very long names
        if len(n) > 30:
            score += len(n) - 30
        # Penalize names with brackets/annotations
        if '(' in n or '[' in n or ',' in n:
            score += 30
        # Tiny penalty for length to break ties
        score += len(n) * 0.01
        return score

    sorted_candidates = sorted(candidates, key=synonym_score)

    # Take top candidates, avoiding similar names
    synonyms = []
    for n in sorted_candidates:
        if len(synonyms) >= max_synonyms:
            break
        n_lower = n.lower()
        # Skip if too similar to primary
        if n_lower in primary_lower or primary_lower in n_lower:
            continue
        # Skip if similar to already added
        if any(n_lower in s.lower() or s.lower() in n_lower for s in synonyms):
            continue
        # Skip code-like or too long
        if len(n) > 40 or synonym_score(n) > 100:
            continue
        synonyms.append(n)

    return synonyms


def _build_basic_result(m, original_query: str) -> Dict[str, Any]:
    """Build a basic result from migrant object (without toxtree)."""
    name = _get_best_name(m, original_query)
    cas = None
    if hasattr(m, 'CAS') and m.CAS:
        cas = m.CAS if isinstance(m.CAS, list) else [m.CAS]

    return {
        "cid": m.cid,
        "name": name,
        "cas": cas,
        "mw": m.M if hasattr(m, 'M') else None,
        "formula": m.formula if hasattr(m, 'formula') else None,
        "smiles": m.smiles if hasattr(m, 'smiles') else None,
        "logP": _extract_scalar(m.logP) if hasattr(m, 'logP') else None,
        "thumbnail_url": f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={m.cid}&t=l",
        "regulatory": {
            "EU": {"authorized": m.EUauthorized if hasattr(m, 'EUauthorized') else None},
            "US": {"authorized": m.USauthorized if hasattr(m, 'USauthorized') else None},
            "CN": {"authorized": m.CNauthorized if hasattr(m, 'CNauthorized') else None},
        },
        "toxtree": None,
    }


def search_pubchem(query: str, max_results: int = 10) -> List[Dict[str, Any]]:
    """
    Search PubChem for substances.

    Uses migrant class for caching and API access.
    """
    loader = get_pubchem_loader()
    if loader is None:
        return []

    try:
        # Try to load by name, CAS, or CID
        result = loader(query, verbose=False)
        if result is not None and hasattr(result, 'cid') and result.cid:
            # Handle name attribute which can be a list
            name = result.name if hasattr(result, 'name') else query
            if isinstance(name, list):
                query_lower = query.lower().strip()

                # First: exact match with user's query (case-insensitive)
                for n in name:
                    if n.lower() == query_lower:
                        name = n
                        break
                else:
                    # Second: user's query contained in a name
                    for n in name:
                        if query_lower in n.lower():
                            name = n
                            break
                    else:
                        # Third: find the best name by scoring
                        def name_score(n):
                            # Lower score = better
                            score = 0
                            if n.startswith(('DTXSID', 'SCHEMBL', 'CHEBI', 'NSC', 'EINECS', 'EC ', 'UNII-', 'AKOS')):
                                score += 100  # Code-like names
                            if any(c.isdigit() for c in n[:3]):
                                score += 50  # Starts with digits (CAS-like)
                            if len(n) > 30:
                                score += len(n) - 30  # Penalize long names
                            if '(' in n or '[' in n:
                                score += 20  # Penalize names with annotations
                            if n.isupper() and len(n) > 5:
                                score += 30  # Penalize all-caps codes
                            return score

                        sorted_names = sorted(name, key=name_score)
                        name = sorted_names[0] if sorted_names else query

            # Handle logP which can be an array
            logP = None
            if hasattr(result, 'logP') and result.logP is not None:
                logP_val = result.logP
                if hasattr(logP_val, '__iter__') and not isinstance(logP_val, str):
                    logP = float(logP_val[0]) if len(logP_val) > 0 else None
                else:
                    logP = float(logP_val)

            return [{
                "cid": result.cid,
                "name": name,
                "cas": result.CAS[0] if hasattr(result, 'CAS') and result.CAS else None,
                "iupac_name": None,  # Not directly available
                "mw": result.M if hasattr(result, 'M') else None,
                "formula": result.formula if hasattr(result, 'formula') else None,
                "smiles": result.smiles if hasattr(result, 'smiles') else None,
                "logP": logP,
                "inchi": result.InChi if hasattr(result, 'InChi') else None,
                "inchikey": result.InChiKey if hasattr(result, 'InChiKey') else None,
                "thumbnail_url": f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={result.cid}&t=l",
            }]
    except Exception as e:
        # Log error but don't fail
        import traceback
        traceback.print_exc()

    return []


def get_substance_by_cid(cid: int) -> Optional[Dict[str, Any]]:
    """Get substance details by PubChem CID using get_compounds with rate limiting."""
    try:
        from patankar.loadpubchem import pubchem_query_with_retry, get_compounds, PubChemQueryError

        # Use rate-limited query for direct CID lookup
        try:
            results = pubchem_query_with_retry(get_compounds, cid, 'cid')
        except PubChemQueryError:
            return None
        if not results:
            return None

        compound = results[0]

        # Get the best name from synonyms
        name = f"CID {cid}"
        if hasattr(compound, 'synonyms') and compound.synonyms:
            synonyms = compound.synonyms
            # Find the best name: prefer short, simple names without codes
            def name_score(n):
                score = 0
                if n.startswith(('DTXSID', 'SCHEMBL', 'CHEBI', 'NSC', 'EINECS', 'EC ', 'UNII-', 'AKOS', 'RefChem')):
                    score += 100
                if any(c.isdigit() for c in n[:3]):
                    score += 50
                if len(n) > 30:
                    score += len(n) - 30
                if '(' in n or '[' in n:
                    score += 20
                if n.isupper() and len(n) > 5:
                    score += 30
                return score
            sorted_names = sorted(synonyms, key=name_score)
            name = sorted_names[0] if sorted_names else f"CID {cid}"

        # Get logP (xlogp property)
        logP = None
        if hasattr(compound, 'xlogp') and compound.xlogp is not None:
            logP = float(compound.xlogp)

        return {
            "cid": compound.cid,
            "name": name,
            "cas": None,  # Not directly available from get_compounds
            "iupac_name": compound.iupac_name if hasattr(compound, 'iupac_name') else None,
            "mw": compound.molecular_weight if hasattr(compound, 'molecular_weight') else None,
            "formula": compound.molecular_formula if hasattr(compound, 'molecular_formula') else None,
            "smiles": compound.canonical_smiles if hasattr(compound, 'canonical_smiles') else None,
            "logP": logP,
            "inchi": compound.inchi if hasattr(compound, 'inchi') else None,
            "inchikey": compound.inchikey if hasattr(compound, 'inchikey') else None,
            "thumbnail_url": f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={cid}&t=l",
        }
    except ImportError:
        pass
    except Exception:
        pass

    return None


# ========== COMMON SUBSTANCES DATABASE ==========

COMMON_SUBSTANCES = {
    "irganox_1076": {
        "id": "irganox_1076",
        "name": "Irganox 1076",
        "cas": "2082-79-3",
        "cid": 40467,
        "mw": 530.87,
        "logP": 12.5,
        "category": "antioxidant",
        "description": "Primary phenolic antioxidant, widely used in polyolefins",
        "typical_C0": {"LDPE": 1000, "HDPE": 800, "PP": 1000},
        "SML": 6.0,
    },
    "irganox_1010": {
        "id": "irganox_1010",
        "name": "Irganox 1010",
        "cas": "6683-19-8",
        "cid": 71777,
        "mw": 1177.63,
        "logP": 23.4,
        "category": "antioxidant",
        "description": "High MW phenolic antioxidant for long-term protection",
        "typical_C0": {"LDPE": 500, "HDPE": 500, "PP": 500},
        "SML": 6.0,
    },
    "bht": {
        "id": "bht",
        "name": "BHT",
        "cas": "128-37-0",
        "cid": 31404,
        "mw": 220.35,
        "logP": 5.1,
        "category": "antioxidant",
        "description": "Butylated hydroxytoluene - Low MW antioxidant, fast migrating",
        "typical_C0": {"LDPE": 100, "HDPE": 100, "PP": 100},
        "SML": 3.0,
    },
    "uvitex_ob": {
        "id": "uvitex_ob",
        "name": "Uvitex OB",
        "cas": "7128-64-5",
        "cid": 23581,
        "mw": 430.56,
        "logP": 8.4,
        "category": "optical_brightener",
        "description": "Optical brightener for plastics",
        "typical_C0": {"LDPE": 50, "HDPE": 50, "PP": 50, "PET": 50},
        "SML": None,
    },
    "deha": {
        "id": "deha",
        "name": "DEHA",
        "cas": "103-23-1",
        "cid": 7641,
        "mw": 370.57,
        "logP": 8.1,
        "category": "plasticizer",
        "description": "Di(2-ethylhexyl) adipate - Plasticizer for PVC and other polymers",
        "typical_C0": {"PVC": 30000},
        "SML": 18.0,
    },
    "dehp": {
        "id": "dehp",
        "name": "DEHP",
        "cas": "117-81-7",
        "cid": 8343,
        "mw": 390.56,
        "logP": 7.6,
        "category": "plasticizer",
        "description": "Di(2-ethylhexyl) phthalate - Phthalate plasticizer (restricted in food contact)",
        "typical_C0": {"PVC": 40000},
        "SML": 1.5,
    },
    "tinuvin_326": {
        "id": "tinuvin_326",
        "name": "Tinuvin 326",
        "cas": "3896-11-5",
        "cid": 19600,
        "mw": 315.80,
        "logP": 6.5,
        "category": "uv_stabilizer",
        "description": "Benzotriazole UV absorber",
        "typical_C0": {"LDPE": 500, "PP": 500, "PET": 200},
        "SML": None,
    },
    "slip_agent_erucamide": {
        "id": "slip_agent_erucamide",
        "name": "Erucamide",
        "cas": "112-84-5",
        "cid": 5365369,
        "mw": 337.58,
        "logP": 8.2,
        "category": "slip_agent",
        "description": "Fatty acid amide slip agent",
        "typical_C0": {"LDPE": 500, "PP": 500},
        "SML": None,
    },
}


# ========== SML DATABASE ==========

SML_DATABASE = {
    "2082-79-3": {"name": "Irganox 1076", "SML": 6.0, "group_SML": None},
    "6683-19-8": {"name": "Irganox 1010", "SML": 6.0, "group_SML": None},
    "128-37-0": {"name": "BHT", "SML": 3.0, "group_SML": None},
    "103-23-1": {"name": "DEHA", "SML": 18.0, "group_SML": None},
    "117-81-7": {"name": "DEHP", "SML": 1.5, "group_SML": "phthalates"},
    "84-74-2": {"name": "DBP", "SML": 0.3, "group_SML": "phthalates"},
    "85-68-7": {"name": "BBP", "SML": 30.0, "group_SML": "phthalates"},
    "131-11-3": {"name": "DMP", "SML": 60.0, "group_SML": None},
}

GROUP_SML = {
    "phthalates": 60.0,  # Combined limit for phthalates group
}


# ========== API ENDPOINTS ==========

@router.get("/search")
async def search_substances(
    q: str = Query(..., min_length=2, description="Search query"),
    source: str = Query("all", description="Source: all, local, pubchem"),
    max_results: int = Query(10, ge=1, le=50),
):
    """
    Search for substances by name, CAS, or CID.

    Returns matching substances from local database and/or PubChem.
    """
    results = []

    # Search local database
    if source in ["all", "local"]:
        q_lower = q.lower()
        for sub_id, sub in COMMON_SUBSTANCES.items():
            if (q_lower in sub["name"].lower() or
                q_lower in sub_id.lower() or
                (sub.get("cas") and q_lower in sub["cas"])):
                results.append({
                    "source": "local",
                    "id": sub_id,
                    **sub,
                })

    # Search PubChem if enabled and we have room for more results
    if source in ["all", "pubchem"] and len(results) < max_results:
        # Collect CAS numbers from local results to avoid duplicates
        local_cas_set = set()
        local_cid_set = set()
        for r in results:
            if r.get("cas"):
                cas = r["cas"]
                if isinstance(cas, list):
                    local_cas_set.update(cas)
                else:
                    local_cas_set.add(cas)
            if r.get("cid"):
                local_cid_set.add(r["cid"])

        pubchem_results = search_pubchem(q, max_results - len(results))
        for r in pubchem_results:
            # Skip if same CAS already in local results
            r_cas = r.get("cas")
            if r_cas:
                if isinstance(r_cas, list):
                    if any(c in local_cas_set for c in r_cas):
                        continue
                elif r_cas in local_cas_set:
                    continue
            # Skip if same CID already in local results
            if r.get("cid") and r["cid"] in local_cid_set:
                continue

            results.append({
                "source": "pubchem",
                "id": f"pubchem_{r.get('cid', 'unknown')}",
                **r,
            })

    return JSONResponse({
        "success": True,
        "query": q,
        "results": results[:max_results],
        "count": len(results[:max_results]),
    })


@router.get("/common")
async def list_common_substances(
    category: Optional[str] = Query(None, description="Filter by category"),
):
    """List commonly used food contact substances."""
    substances = list(COMMON_SUBSTANCES.values())

    if category:
        substances = [s for s in substances if s.get("category") == category]

    categories = list(set(s.get("category") for s in COMMON_SUBSTANCES.values() if s.get("category")))

    return JSONResponse({
        "success": True,
        "substances": substances,
        "count": len(substances),
        "categories": sorted(categories),
    })


@router.get("/pubchem/{cid}")
async def get_pubchem_substance(cid: int):
    """Get substance details from PubChem by CID."""
    result = get_substance_by_cid(cid)

    if result is None:
        raise HTTPException(status_code=404, detail=f"Substance CID {cid} not found")

    return JSONResponse({
        "success": True,
        "substance": result,
    })


@router.get("/detail/{query}")
async def get_substance_detail(query: str):
    """
    Get comprehensive substance details including toxicological data.

    Uses loadpubchem.migrant and migrantToxtree for full information:
    - PubChem properties (MW, logP, SMILES, etc.)
    - EU/US/CN regulatory authorization status
    - SML values
    - Cramer toxicological classification
    - Toxicological structural alerts
    """
    result = load_substance_with_toxtree(query)

    if result is None:
        # Try basic search
        pubchem_results = search_pubchem(query, 1)
        if pubchem_results:
            return JSONResponse({
                "success": True,
                "substance": pubchem_results[0],
                "full_data": False,
                "message": "Basic data only (migrant class not available)",
            })
        raise HTTPException(
            status_code=404,
            detail=f"Substance '{query}' not found in PubChem. Use a valid CAS number, IUPAC name, brand name (e.g., Irganox 1076), or known acronym (e.g., BHT)."
        )

    return JSONResponse({
        "success": True,
        "substance": result,
        "full_data": True,
        "has_toxtree": result.get("toxtree") is not None,
    })


@router.get("/regulatory/{query}")
async def get_regulatory_status(query: str):
    """
    Get regulatory authorization status for a substance.

    Returns authorization status in:
    - EU (EU 10/2011)
    - US (FDA 21 CFR)
    - CN (GB 9685)
    """
    result = load_substance_with_toxtree(query)

    if result is None:
        raise HTTPException(
            status_code=404,
            detail=f"Substance '{query}' not found in PubChem. Use a valid CAS number, IUPAC name, brand name (e.g., Irganox 1076), or known acronym (e.g., BHT)."
        )

    return JSONResponse({
        "success": True,
        "substance": result.get("name"),
        "cid": result.get("cid"),
        "cas": result.get("cas"),
        "regulatory": result.get("regulatory", {}),
    })


@router.get("/sml/{cas}")
async def get_sml(cas: str):
    """Get SML (Specific Migration Limit) for a substance by CAS number."""
    if cas in SML_DATABASE:
        info = SML_DATABASE[cas]
        group_limit = None
        if info.get("group_SML"):
            group_limit = GROUP_SML.get(info["group_SML"])

        return JSONResponse({
            "success": True,
            "cas": cas,
            "name": info["name"],
            "SML": info["SML"],
            "unit": "mg/kg food",
            "group_SML": info.get("group_SML"),
            "group_limit": group_limit,
            "regulation": "EU 10/2011",
        })

    return JSONResponse({
        "success": False,
        "cas": cas,
        "message": "SML not found in database",
        "SML": None,
    }, status_code=404)


@router.post("/validate")
async def validate_substance_config(substances: List[SubstanceInLayer]):
    """Validate substance configuration for a simulation."""
    errors = []
    warnings = []
    validated = []

    for sub in substances:
        entry = {"substance_id": sub.substance_id, "layer_index": sub.layer_index}

        # Check if substance exists
        if sub.substance_id in COMMON_SUBSTANCES:
            info = COMMON_SUBSTANCES[sub.substance_id]
            entry["name"] = info["name"]
            entry["mw"] = info.get("mw")
            entry["logP"] = info.get("logP")

            # Check SML
            if sub.use_SML_check:
                if sub.SML is None and info.get("SML"):
                    entry["SML"] = info["SML"]
                    warnings.append(f"{info['name']}: Using default SML = {info['SML']} mg/kg")
                elif sub.SML is None:
                    warnings.append(f"{info['name']}: No SML defined, compliance check disabled")
                    entry["SML"] = None
                else:
                    entry["SML"] = sub.SML
        else:
            # Unknown substance - need MW and logP for estimation
            entry["name"] = sub.substance_id
            warnings.append(f"Unknown substance '{sub.substance_id}': provide MW and logP for D estimation")

        # Validate C0
        if sub.C0 <= 0:
            errors.append(f"Invalid C0 for {sub.substance_id}: must be > 0")
        elif sub.C0 > 100000:
            warnings.append(f"Very high C0 for {sub.substance_id}: {sub.C0} mg/kg")

        entry["C0"] = sub.C0
        validated.append(entry)

    return JSONResponse({
        "success": len(errors) == 0,
        "valid": len(errors) == 0,
        "errors": errors,
        "warnings": warnings,
        "substances": validated,
    })


@router.get("/categories")
async def list_substance_categories():
    """List substance categories."""
    categories = {
        "antioxidant": {
            "code": "antioxidant",
            "name": "Antioxidants",
            "icon": "🛡️",
            "description": "Primary and secondary antioxidants for polymer stabilization",
            "examples": ["Irganox 1076", "Irganox 1010", "BHT"],
        },
        "plasticizer": {
            "code": "plasticizer",
            "name": "Plasticizers",
            "icon": "💧",
            "description": "Softening agents for flexible polymers",
            "examples": ["DEHA", "DEHP", "DINP"],
        },
        "uv_stabilizer": {
            "code": "uv_stabilizer",
            "name": "UV Stabilizers",
            "icon": "☀️",
            "description": "Light stabilizers and UV absorbers",
            "examples": ["Tinuvin 326", "Tinuvin 770"],
        },
        "slip_agent": {
            "code": "slip_agent",
            "name": "Slip Agents",
            "icon": "🧴",
            "description": "Surface modifiers for reduced friction",
            "examples": ["Erucamide", "Oleamide"],
        },
        "optical_brightener": {
            "code": "optical_brightener",
            "name": "Optical Brighteners",
            "icon": "✨",
            "description": "Fluorescent whitening agents",
            "examples": ["Uvitex OB"],
        },
        "monomer": {
            "code": "monomer",
            "name": "Residual Monomers",
            "icon": "🔬",
            "description": "Unreacted monomers from polymerization",
            "examples": ["Styrene", "Vinyl chloride", "Acrylonitrile"],
        },
    }

    return JSONResponse({
        "success": True,
        "categories": list(categories.values()),
    })


@router.get("/thumbnail-urls/{cid}")
async def get_substance_thumbnail_urls(cid: int):
    """Get thumbnail URLs for a PubChem substance (metadata only, no image)."""
    return JSONResponse({
        "success": True,
        "cid": cid,
        "thumbnail_url": f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={cid}&t=l",
        "structure_2d_url": f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={cid}&t=s",
        "structure_3d_url": f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG?record_type=3d",
        "local_tiny_url": f"/api/substances/thumbnail/{cid}.png",
        "local_large_url": f"/api/substances/structure/{cid}.png",
    })


@router.get("/debug/toxtree/{query}")
async def debug_toxtree(query: str):
    """
    Debug endpoint for ToxTree data inspection.

    Returns raw migrantToxtree data for troubleshooting.
    """
    import traceback

    result = {
        "query": query,
        "migrant_available": False,
        "toxtree_available": False,
        "migrant_data": None,
        "toxtree_data": None,
        "errors": [],
    }

    # Try basic migrant first
    migrant_cls = get_migrant_class()
    if migrant_cls:
        result["migrant_available"] = True
        try:
            m = migrant_cls(query, verbose=False)
            if m and hasattr(m, 'cid') and m.cid:
                result["migrant_data"] = {
                    "cid": m.cid,
                    "name": m.name if hasattr(m, 'name') else None,
                    "smiles": m.smiles if hasattr(m, 'smiles') else None,
                    "cas": m.CAS if hasattr(m, 'CAS') else None,
                }
        except Exception as e:
            result["errors"].append(f"migrant error: {str(e)}")

    # Try migrantToxtree
    toxtree_cls = get_migrant_toxtree_class()
    if toxtree_cls:
        result["toxtree_available"] = True
        try:
            mt = toxtree_cls(query, verbose=False)
            if mt and hasattr(mt, 'cid') and mt.cid:
                # Basic info
                toxtree_info = {
                    "cid": mt.cid,
                    "name": mt.name if hasattr(mt, 'name') else None,
                    "smiles": mt.smiles if hasattr(mt, 'smiles') else None,
                }

                # Cramer classification
                if hasattr(mt, 'CramerClass'):
                    toxtree_info["cramer_class"] = mt.CramerClass
                    toxtree_info["cramer_value"] = getattr(mt, 'CramerValue', None)
                    toxtree_info["ttc"] = getattr(mt, 'TTC', None)
                    toxtree_info["ttc_units"] = getattr(mt, 'TTCunits', None)
                    toxtree_info["cf_ttc"] = getattr(mt, 'CFTTC', None)

                # Alerts from skin plugin
                toxtree_info["has_alerts"] = getattr(mt, 'has_alerts', False)
                toxtree_info["nalerts"] = getattr(mt, 'nalerts', 0)
                toxtree_info["showalerts"] = getattr(mt, 'showalerts', {})

                # DNA binding check
                try:
                    if hasattr(mt, '_run_toxtree'):
                        dna_result = mt._run_toxtree('dnabinding')
                        dna_alerts = {}
                        for key in ['AlertForMichaelAcceptorIdentified', 'AlertForAcylTransferAgentIdentified',
                                    'AlertForSn1Identified', 'AlertForSn2Identified', 'AlertForSchiffBaseFormationIdentified']:
                            dna_alerts[key] = dna_result.get(key, 'N/A')
                        toxtree_info["dna_binding_alerts"] = dna_alerts
                except Exception as e:
                    toxtree_info["dna_binding_error"] = str(e)

                result["toxtree_data"] = toxtree_info
        except Exception as e:
            result["errors"].append(f"toxtree error: {str(e)}\n{traceback.format_exc()}")

    return JSONResponse(result)


# ========== OFFLINE STRUCTURE IMAGES ==========

# Path to the PubChem cache directory
PUBCHEM_CACHE_DIR = Path(__file__).parent.parent.parent.parent / "patankar" / "cache.PubChem"
# Large images (800x800) - managed by loadpubchem.py
PUBCHEM_THUMBS_DIR = PUBCHEM_CACHE_DIR / "thumbs"
# Tiny thumbnails (~150px) - for studio UI, separate from loadpubchem
PUBCHEM_TINY_DIR = PUBCHEM_CACHE_DIR / "thumbs_tiny"


@router.get("/structure/available")
async def list_cached_structures():
    """
    List all CIDs with cached structure images (both large and tiny).

    Useful for offline mode to know which substances have local images.
    """
    # Large images (800x800) from loadpubchem
    large_cids = []
    if PUBCHEM_THUMBS_DIR.exists():
        for f in PUBCHEM_THUMBS_DIR.glob("*.png"):
            try:
                cid = int(f.stem)
                large_cids.append(cid)
            except ValueError:
                pass

    # Tiny thumbnails (~150px) for studio UI
    tiny_cids = []
    if PUBCHEM_TINY_DIR.exists():
        for f in PUBCHEM_TINY_DIR.glob("*.png"):
            try:
                cid = int(f.stem)
                tiny_cids.append(cid)
            except ValueError:
                pass

    return JSONResponse({
        "success": True,
        "large": {
            "count": len(large_cids),
            "cached_cids": sorted(large_cids),
            "cache_dir": str(PUBCHEM_THUMBS_DIR)
        },
        "tiny": {
            "count": len(tiny_cids),
            "cached_cids": sorted(tiny_cids),
            "cache_dir": str(PUBCHEM_TINY_DIR)
        }
    })


@router.get("/structure/{cid}.png", summary="Get structure image with auto-caching")
async def get_structure_image_with_cache(cid: int, auto_cache: bool = True):
    """
    Serve structure image from local cache, with optional auto-download from PubChem.

    This endpoint:
    1. First checks local cache (patankar/cache.PubChem/thumbs/)
    2. If not cached and auto_cache=True, downloads from PubChem and caches locally
    3. Returns 404 if not cached and auto_cache=False

    This enables offline operation by caching images as they're accessed.

    Parameters:
    - cid: PubChem Compound ID
    - auto_cache: If True, download and cache images not found locally (default: True)
    """
    import httpx

    image_path = PUBCHEM_THUMBS_DIR / f"{cid}.png"

    # 1. Check local cache first
    if image_path.exists():
        return FileResponse(
            path=str(image_path),
            media_type="image/png",
            headers={
                "Cache-Control": "public, max-age=31536000",  # 1 year cache
                "X-Source": "local-cache"
            }
        )

    # 2. If auto_cache enabled, try to download and cache
    if auto_cache:
        # PubChem image URL - use larger size for better quality
        pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={cid}&t=l"

        try:
            async with httpx.AsyncClient(timeout=10.0) as client:
                response = await client.get(pubchem_url)
                response.raise_for_status()

                # Ensure thumbs directory exists
                PUBCHEM_THUMBS_DIR.mkdir(parents=True, exist_ok=True)

                # Save to cache
                image_path.write_bytes(response.content)

                return Response(
                    content=response.content,
                    media_type="image/png",
                    headers={
                        "Cache-Control": "public, max-age=31536000",
                        "X-Source": "auto-cached-from-pubchem"
                    }
                )
        except httpx.HTTPError as e:
            # Log error but don't fail - return 404 with fallback URL
            print(f"[Substances] Failed to auto-cache CID {cid}: {e}")
        except Exception as e:
            print(f"[Substances] Unexpected error caching CID {cid}: {e}")

    # 3. Not cached and couldn't download
    raise HTTPException(
        status_code=404,
        detail={
            "error": "Structure image not cached",
            "cid": cid,
            "hint": "Use online mode to cache this substance first",
            "fallback_url": f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={cid}&t=l"
        }
    )


@router.post("/structure/cache-batch")
async def cache_structure_images_batch(cids: List[int]):
    """
    Cache structure images for multiple CIDs in batch.

    Useful for pre-caching images for offline use.

    Parameters:
    - cids: List of PubChem CIDs to cache (max 50)

    Returns:
    - success: list of CIDs successfully cached
    - already_cached: list of CIDs already in cache
    - failed: list of CIDs that failed to cache
    """
    import httpx

    if len(cids) > 50:
        raise HTTPException(
            status_code=400,
            detail="Maximum 50 CIDs per batch request"
        )

    results = {
        "success": [],
        "already_cached": [],
        "failed": []
    }

    # Ensure thumbs directory exists
    PUBCHEM_THUMBS_DIR.mkdir(parents=True, exist_ok=True)

    async with httpx.AsyncClient(timeout=10.0) as client:
        for cid in cids:
            image_path = PUBCHEM_THUMBS_DIR / f"{cid}.png"

            # Skip if already cached
            if image_path.exists():
                results["already_cached"].append(cid)
                continue

            # Try to download
            try:
                pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={cid}&t=l"
                response = await client.get(pubchem_url)
                response.raise_for_status()

                image_path.write_bytes(response.content)
                results["success"].append(cid)

                # Small delay to be nice to PubChem
                import asyncio
                await asyncio.sleep(0.2)

            except Exception as e:
                results["failed"].append({"cid": cid, "error": str(e)})

    return JSONResponse({
        "success": True,
        "cached": len(results["success"]),
        "already_cached": len(results["already_cached"]),
        "failed": len(results["failed"]),
        "details": results
    })


# ========== TINY THUMBNAILS (for studio UI) ==========

@router.get("/thumbnail/{cid}.png", summary="Get tiny thumbnail with auto-caching")
async def get_tiny_thumbnail_with_cache(cid: int, auto_cache: bool = True):
    """
    Serve tiny thumbnail (~150px) from local cache, with optional auto-download from PubChem.

    This endpoint serves the SMALL thumbnails used in the studio UI (not the 800x800 images
    from loadpubchem). These are cached separately in cache.PubChem/thumbs_tiny/.

    This endpoint:
    1. First checks local tiny cache (patankar/cache.PubChem/thumbs_tiny/)
    2. If not cached and auto_cache=True, downloads small thumbnail from PubChem
    3. Returns 404 if not cached and auto_cache=False

    Parameters:
    - cid: PubChem Compound ID
    - auto_cache: If True, download and cache images not found locally (default: True)
    """
    import httpx

    image_path = PUBCHEM_TINY_DIR / f"{cid}.png"

    # 1. Check local tiny cache first
    if image_path.exists():
        return FileResponse(
            path=str(image_path),
            media_type="image/png",
            headers={
                "Cache-Control": "public, max-age=31536000",  # 1 year cache
                "X-Source": "local-tiny-cache"
            }
        )

    # 2. If auto_cache enabled, try to download small thumbnail
    if auto_cache:
        # PubChem small thumbnail URL (t=s for small, ~150px)
        pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={cid}&t=s"

        try:
            async with httpx.AsyncClient(timeout=10.0) as client:
                response = await client.get(pubchem_url)
                response.raise_for_status()

                # Ensure tiny thumbs directory exists
                PUBCHEM_TINY_DIR.mkdir(parents=True, exist_ok=True)

                # Save to tiny cache
                image_path.write_bytes(response.content)

                return Response(
                    content=response.content,
                    media_type="image/png",
                    headers={
                        "Cache-Control": "public, max-age=31536000",
                        "X-Source": "auto-cached-tiny-from-pubchem"
                    }
                )
        except httpx.HTTPError as e:
            print(f"[Substances] Failed to auto-cache tiny thumbnail CID {cid}: {e}")
        except Exception as e:
            print(f"[Substances] Unexpected error caching tiny CID {cid}: {e}")

    # 3. Not cached and couldn't download
    raise HTTPException(
        status_code=404,
        detail={
            "error": "Tiny thumbnail not cached",
            "cid": cid,
            "hint": "Use online mode to cache this substance first",
            "fallback_url": f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={cid}&t=s"
        }
    )


@router.post("/thumbnail/cache-batch")
async def cache_tiny_thumbnails_batch(cids: List[int]):
    """
    Cache tiny thumbnails for multiple CIDs in batch.

    Useful for pre-caching small thumbnails for offline use in the studio UI.

    Parameters:
    - cids: List of PubChem CIDs to cache (max 50)

    Returns:
    - success: list of CIDs successfully cached
    - already_cached: list of CIDs already in cache
    - failed: list of CIDs that failed to cache
    """
    import httpx

    if len(cids) > 50:
        raise HTTPException(
            status_code=400,
            detail="Maximum 50 CIDs per batch request"
        )

    results = {
        "success": [],
        "already_cached": [],
        "failed": []
    }

    # Ensure tiny thumbs directory exists
    PUBCHEM_TINY_DIR.mkdir(parents=True, exist_ok=True)

    async with httpx.AsyncClient(timeout=10.0) as client:
        for cid in cids:
            image_path = PUBCHEM_TINY_DIR / f"{cid}.png"

            # Skip if already cached
            if image_path.exists():
                results["already_cached"].append(cid)
                continue

            # Try to download small thumbnail
            try:
                pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/image/imgsrv.fcgi?cid={cid}&t=s"
                response = await client.get(pubchem_url)
                response.raise_for_status()

                image_path.write_bytes(response.content)
                results["success"].append(cid)

                # Small delay to be nice to PubChem
                import asyncio
                await asyncio.sleep(0.2)

            except Exception as e:
                results["failed"].append({"cid": cid, "error": str(e)})

    return JSONResponse({
        "success": True,
        "cached": len(results["success"]),
        "already_cached": len(results["already_cached"]),
        "failed": len(results["failed"]),
        "details": results
    })
