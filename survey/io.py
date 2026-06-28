"""
survey/io.py — Input/Output Functions for Survey
================================================

Handles scenario loading from YAML and result dumping to NPZ/JSON.

@project: SFPPy — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

import json
import warnings
from pathlib import Path
from typing import Dict, Any, Union, List, Optional
import logging

logger = logging.getLogger(__name__)

import numpy as np
import yaml

from survey.models import (
    LayerSpec,
    PackagingSpec,
    SubstanceSpec,
    PriorSpec,
    SurveyConfig,
    ComponentSpec,
    MultiComponentJob,
)


# =============================================================================
# Surrogate CAS Mapping for Non-Pure Compounds
# =============================================================================
# When a CAS lookup fails (polymer, mixture, oil), use these surrogates.
# For oils: use TAG (triglyceride) as primary - better MW representation.
# Format: "original_CAS": "surrogate_CAS"  or  "original_CAS": ("TAG_CAS", "FA_CAS")

SURROGATE_CAS = {
    # Polymers → Monomers
    "9002-88-4":  "74-85-1",      # PE (Polyethylene) → Ethylene
    "65447-77-0": "2226-96-2",    # Tinuvin 622 LD → 4-Hydroxy-TEMPO
    "9005-32-7":  "6814-36-4",    # Alginic acid → D-Mannuronic acid
    "9036-19-5":  "140-66-9",     # PEG octylphenyl ether → 4-tert-Octylphenol

    # Flame retardants → Base structure
    "68631-49-2": "101-84-8",     # Hexabromodiphenyl ether → Diphenyl ether

    # Wood tar/Creosote
    "8021-39-4":  "90-05-1",      # Wood creosote → Guaiacol

    # Waxes
    "8002-53-7":  "506-48-9",     # Montan wax → Montanic acid

    # Vegetable oils → (Triglyceride, Fatty acid)
    # Primary: use TAG for better MW; fallback: FA
    "8002-75-3":  ("555-44-2", "57-10-3"),    # Palm oil → Tripalmitin, Palmitic acid
    "8001-30-7":  ("537-40-6", "60-33-3"),    # Corn oil → Trilinolein, Linoleic acid
    "8001-79-4":  ("2540-54-7", "141-22-0"),  # Castor oil → Triricinolein, Ricinoleic acid
    "8002-13-9":  ("122-32-7", "112-80-1"),   # Rapeseed oil → Triolein, Oleic acid
    "8012-95-1":  "544-76-3",                  # Mineral oil → Hexadecane
    "84836-98-6": ("538-24-9", "143-07-7"),   # Hydrogenated coconut → Trilaurin, Lauric acid
    "8001-78-3":  ("555-43-1", "57-11-4"),    # Hydrogenated castor → Tristearin, Stearic acid

    # Surfactants
    "26658-19-5": "57-11-4",      # Sorbitan tristearate → Stearic acid
    "68412-54-4": "104-40-5",     # Nonylphenol ethoxylated → 4-Nonylphenol

    # Phthalates
    "71888-89-6": "117-81-7",     # Branched phthalate → DEHP

    # Carbohydrates
    "68412-29-3": "50-99-7",      # Hydrolyzed starch → D-Glucose

    # Resins
    "65997-06-0": "79-54-9",      # Hydrogenated rosin → Dihydroabietic acid

    # Fatty acid derivatives
    "61791-31-9": "120-40-1",     # Cocamide DEA → Lauramide DEA

    # Acrylates
    "12542-30-2": "1330-61-6",    # Dicyclopentenyl acrylate → Isodecyl acrylate

    # Colorants
    "1325-82-2":  "548-62-9",     # Pigment Violet 3 → Crystal violet

    # ------------------------------------------------------------------
    # UVCB / mixture representatives — dbs_20260618 (3a) paper-overpack
    # families. Each is a NEUTRAL, single-CID, conservative representative
    # (no salts/ions); oils/surfactants→fatty acid, polymer→monomer,
    # rosin→abietic, halogenated→retains halogen (Rev1 ladder rules).
    # Approved OV 2026-06-27.
    # ------------------------------------------------------------------
    "1330-80-9":  "112-80-1",                 # PG monooleate → Oleic acid
    "1336-36-3":  "2051-60-7",                # PCBs → 2-Chlorobiphenyl (retains Cl)
    "1401-55-4":  "149-91-7",                 # Tannic acid → Gallic acid
    "26266-58-0": "112-80-1",                 # Sorbitan trioleate → Oleic acid
    "28804-88-8": "581-42-0",                 # Dimethylnaphthalene → 2,6-Dimethylnaphthalene
    "54276-35-6": "80-62-6",                  # Sulfopropyl methacrylate → Methyl methacrylate
    "61791-12-6": "141-22-0",                 # PEG-30 castor oil → Ricinoleic acid
    "63148-62-9": "556-67-2",                 # PDMS → Octamethylcyclotetrasiloxane (D4)
    "67762-27-0": "36653-82-4",               # Cetearyl alcohol → 1-Hexadecanol (cetyl)
    "68441-17-8": "57-11-4",                  # Oxidised PE wax → Stearic acid
    "68610-06-0": "88-18-6",                  # Isobutylenated phenol → 2-tert-Butylphenol
    "73138-82-6": "514-10-3",                 # Resin acids → Abietic acid
    "8001-22-7":  ("537-40-6", "60-33-3"),    # Soybean oil → Trilinolein, Linoleic acid
    "8006-44-8":  "630-04-6",                 # Candelilla wax → Hentriacontane (C31)
    "8012-89-3":  "57-10-3",                  # Beeswax → Palmitic acid
    "85535-85-9": "58-89-9",                  # MCCPs → Lindane (retains Cl)

    # Fluorescent whitening agents (sulfonated stilbenes, M~1300, logP missing
    # → non-finite partition → solver abort). Conservative neutral stilbene-core
    # representative (trans-stilbene: smallest/most-mobile, valid logP). OV 2026-06-27.
    "41098-56-0": "103-30-0",                 # Hexasulfo FWA → trans-Stilbene
    "68971-49-3": "103-30-0",                 # Fluorescent Brightener 264 → trans-Stilbene
}


def _migrant_usable(m) -> bool:
    """
    True if a resolved migrant carries the descriptors the survey needs to
    build a finite initial state.

    Some substances resolve to a valid PubChem CID yet lack a usable ``logP``
    (None / non-finite) — e.g. high-MW sulfonated fluorescent whitening agents
    (M~1300). The food/polymer partition step then yields NaN and the solver
    aborts ("initial state y0 must be finite"). Such a substance must be
    treated like an unresolved one so the SURROGATE_CAS ladder can supply a
    computable representative.
    """
    logP = getattr(m, "logP", None)
    if logP is None:
        return False
    try:
        arr = np.ravel(np.asarray(logP, dtype=float))
    except (TypeError, ValueError):
        return False
    return arr.size > 0 and bool(np.all(np.isfinite(arr)))


def _get_surrogate_cas(cas: str) -> Optional[List[str]]:
    """
    Get surrogate CAS number(s) for a failed lookup.

    Parameters
    ----------
    cas : str
        Original CAS number that failed.

    Returns
    -------
    list of str or None
        List of surrogate CAS numbers to try, or None if no surrogate exists.
        For oils, returns [TAG_CAS, FA_CAS] to try both.
    """
    if cas not in SURROGATE_CAS:
        return None

    surrogate = SURROGATE_CAS[cas]
    if isinstance(surrogate, tuple):
        return list(surrogate)  # TAG and FA for oils
    return [surrogate]


def load_yaml(path: Union[str, Path]) -> Dict[str, Any]:
    """
    Load YAML file.

    Parameters
    ----------
    path : str or Path
        Path to YAML file.

    Returns
    -------
    dict
        Parsed YAML content.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    with open(path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def parse_layer_specs(
    physics: Dict[str, Any],
    interface_contact_temperature_degC: Optional[float] = None,
) -> List[LayerSpec]:
    """
    Parse layer specifications from physics config.

    Supports both monolayer and multilayer formats.

    Temperature priority (most specific first) — Patch P2 (2026-05-12):
      1. ``physics.{mono,multi}.layers[i].temperature_degC``  (per-layer)
      2. ``physics.{mono,multi}.temperature_degC``            (parent-level)
      3. ``physics.interface.contact_temperature_degC``       (passed in by caller)
      4. 40.0                                                  (last-resort default)

    Before P2, the multilayer branch read only the parent-level value
    and applied it to every layer; the per-layer field was silently
    ignored (fixed in the v3 survey-engine consolidation).

    Parameters
    ----------
    physics : dict
        Physics section from scenario YAML.
    interface_contact_temperature_degC : float, optional
        Fallback temperature when neither the per-layer nor the
        parent-level field is present.  Threaded in from
        ``parse_packaging_spec``; preserves backward compat with
        scenarios that omit per-layer T.

    Returns
    -------
    List[LayerSpec]
        List of layer specifications.
    """
    iface_T = (
        float(interface_contact_temperature_degC)
        if interface_contact_temperature_degC is not None
        else 40.0
    )

    # Check for monolayer format
    if 'monolayer' in physics:
        mono = physics['monolayer']
        return [LayerSpec(
            polymer=mono.get('polymer', 'LDPE'),
            thickness_m=float(mono.get('thickness_m', 100e-6)),
            C0=float(mono.get('C0', 0.0)),
            temperature_degC=float(mono.get('temperature_degC', iface_T)),
        )]

    # Check for multilayer format
    if 'multilayer' in physics:
        multi = physics['multilayer']
        layers_cfg = multi.get('layers', [])
        multi_T = float(multi.get('temperature_degC', iface_T))
        return [
            LayerSpec(
                polymer=layer.get('polymer', 'LDPE'),
                thickness_m=float(layer.get('thickness_m', 100e-6)),
                C0=float(layer.get('C0', 0.0)),
                temperature_degC=float(layer.get('temperature_degC', multi_T)),
            )
            for layer in layers_cfg
        ]

    raise ValueError("Physics config must have 'monolayer' or 'multilayer' section")


def parse_packaging_spec(physics: Dict[str, Any]) -> PackagingSpec:
    """
    Parse full packaging specification from physics config.

    Parameters
    ----------
    physics : dict
        Physics section from scenario YAML.

    Returns
    -------
    PackagingSpec
        Packaging specification.
    """
    iface = physics.get('interface', {})
    iface_T = float(iface.get('contact_temperature_degC', 40.0))

    # Patch P2 (2026-05-12): thread the contact temperature into
    # parse_layer_specs as the per-layer fallback when neither
    # `layers[i].temperature_degC` nor `multilayer.temperature_degC`
    # is provided.
    layers = parse_layer_specs(physics, interface_contact_temperature_degC=iface_T)

    return PackagingSpec(
        layers=layers,
        h_m_s=float(iface.get('h_m_s', 1e-7)),
        surface_area_m2=float(iface.get('surface_area_m2', 0.06)),
        food_volume_m3=float(iface.get('food_volume_m3', 0.001)),
        contact_temperature_degC=iface_T,
        cf0=float(iface.get('cf0', 0.0)),
        food_simulant=iface.get('food_simulant', 'oliveoil'),
    )


def parse_prior_spec(prior_cfg: Dict[str, Any], name: str = "prior") -> PriorSpec:
    """
    Parse prior specification from config.

    Parameters
    ----------
    prior_cfg : dict
        Prior configuration with 'triangular' and 'grid' sections.
    name : str
        Name for the prior.

    Returns
    -------
    PriorSpec
        Prior specification.
    """
    tri = prior_cfg.get('triangular', {})
    grid = prior_cfg.get('grid', {})

    return PriorSpec(
        min_val=float(tri.get('min', 0.0)),
        mode=float(tri.get('mode', 1.0)),
        max_val=float(tri.get('max', 10.0)),
        n_low=int(grid.get('nlow', 15)),
        n_high=int(grid.get('nhigh', 15)),
        spacing_low=str(grid.get('spacing_low', 'linear')),
        spacing_high=str(grid.get('spacing_high', 'linear')),
        name=name,
    )


def _try_migrant_lookup(identifier: str) -> tuple:
    """
    Try to resolve a substance via migrant(), checking for non-pure compounds.

    Parameters
    ----------
    identifier : str
        CAS number or substance name to look up.

    Returns
    -------
    tuple
        (migrant_obj, molar_mass) if successful.

    Raises
    ------
    ValueError
        If lookup fails or returns non-pure compound.
    """
    from patankar.loadpubchem import migrant

    m = migrant(name=identifier)
    M = float(m.M)

    # Check if resolved to a valid molar mass
    if M <= 0 or M > 50000:
        raise ValueError(f"Invalid molar mass {M} for '{identifier}'")

    return m, M


def parse_substances(
    family_cfg: Dict[str, Any],
    *,
    skip_unresolvable: bool = True,
    data_gaps: Optional[List[Dict[str, Any]]] = None,
) -> List[SubstanceSpec]:
    """
    Parse substance specifications from family config.

    Supports:
    1. Direct masses: family.masses_g_mol = [80, 120, ...]
    2. Named substances: family.substances = [{name: "...", cas: "..."}, ...]

    Resolution ladder (per substance): original CAS → name → surrogate CAS
    (SURROGATE_CAS, the curated "representative conservative substance" map
    for polymers, oils, waxes, mixtures).

    Exchangeability semantics: a family's ``substances`` list is an
    *uncertainty set* — each member is simulated individually and combined
    by a mixture rule downstream (R5b). A member that cannot be resolved
    (a UVCB / mixture with no single PubChem CID and no ladder entry, e.g.
    soybean oil, beeswax, PCBs, PDMS) therefore must NOT silently disappear
    nor abort the whole set. With ``skip_unresolvable=True`` (default) such a
    member is reported LOUDLY (``logger.error`` + ``warnings.warn``), recorded
    as a flagged data-gap, and skipped — the resolvable members still run.
    A family in which EVERY member is unresolvable still raises (a true,
    surfaced failure rather than an empty silent result).

    Parameters
    ----------
    family_cfg : dict
        Family section from scenario YAML.
    skip_unresolvable : bool, keyword-only, default True
        If True, accept unresolvable members as flagged data-gaps and
        continue (the curated-ladder + accept-failures policy). If False,
        restore the strict contract: raise on the first unresolvable member.
    data_gaps : list, keyword-only, optional
        If provided, each skipped member is appended as
        ``{"cas", "name", "reason"}`` for downstream persistence/reporting.

    Returns
    -------
    List[SubstanceSpec]
        List of substance specifications (resolvable members only when
        ``skip_unresolvable`` is True).
    """
    # Direct masses
    if 'masses_g_mol' in family_cfg:
        masses = family_cfg['masses_g_mol']
        return [
            SubstanceSpec.from_mass(float(m), idx=i)
            for i, m in enumerate(masses)
        ]

    # Named substances
    if 'substances' in family_cfg:
        substances = []
        for i, spec in enumerate(family_cfg['substances']):
            name = spec.get('name')
            cas = spec.get('cas')

            if name is None and cas is None:
                raise ValueError(f"Substance spec must have 'name' or 'cas': {spec}")

            # Resolution strategy:
            # 1. Try original CAS (if provided)
            # 2. Try name as fallback (if CAS failed and name provided)
            # 3. Try surrogate CAS numbers (if original CAS in SURROGATE_CAS)
            # 4. Raise error if all attempts fail

            m = None
            M = None
            used_surrogate = None
            original_error = None

            # Step 1: Try original CAS or name
            try:
                m, M = _try_migrant_lookup(cas or name)
            except Exception as e:
                original_error = e

                # Step 2: If CAS failed and we have a name, try name
                if cas and name:
                    try:
                        m, M = _try_migrant_lookup(name)
                        logger.debug(f"CAS '{cas}' failed, resolved via name '{name}'")
                    except Exception:
                        pass  # Continue to surrogate lookup

                # Step 3: Try surrogate CAS numbers
                if m is None and cas:
                    surrogates = _get_surrogate_cas(cas)
                    if surrogates:
                        for surrogate_cas in surrogates:
                            try:
                                m, M = _try_migrant_lookup(surrogate_cas)
                                used_surrogate = surrogate_cas
                                logger.warning(
                                    f"Using surrogate CAS '{surrogate_cas}' for "
                                    f"'{cas}' ({name or 'unnamed'})"
                                )
                                break
                            except Exception:
                                continue

            # Step 3.5: A substance can RESOLVE yet lack the descriptors needed
            # to simulate (e.g. logP=None for high-MW sulfonated FWAs → NaN
            # partition → solver abort). If — and ONLY if — a usable surrogate
            # is defined for it, redirect to that surrogate. Otherwise KEEP the
            # originally-resolved migrant unchanged: many logP-less substances
            # (metals, salts: As, Ba, Zn stearate…) are still handled by the
            # engine, and were computed by the pre-gate code — this gate must
            # NOT drop them (no regression). It is a narrow surrogate-redirect,
            # never a discard.
            if m is not None and used_surrogate is None and not _migrant_usable(m):
                surrogates = _get_surrogate_cas(cas) if cas else None
                if surrogates:
                    for surrogate_cas in surrogates:
                        try:
                            m_s, M_s = _try_migrant_lookup(surrogate_cas)
                        except Exception:
                            continue
                        if _migrant_usable(m_s):
                            m, M, used_surrogate = m_s, M_s, surrogate_cas
                            logger.warning(
                                f"Substance '{cas}' ({name or 'unnamed'}) resolved but "
                                f"lacks usable descriptors (logP); using surrogate CAS "
                                f"'{surrogate_cas}'."
                            )
                            break

            # Step 4: If all attempts failed — either raise (strict) or
            # accept the member as a LOUD, recorded data-gap (default).
            if m is None:
                surrogates = _get_surrogate_cas(cas) if cas else None
                if surrogates:
                    reason = (
                        f"Failed to resolve substance '{cas}' (CAS) or '{name}' (name), "
                        f"and surrogates {surrogates} also failed: {original_error}"
                    )
                else:
                    reason = f"Failed to resolve substance '{cas or name}': {original_error}"
                if not skip_unresolvable:
                    raise ValueError(reason)
                # Accept-failure policy: the member is a UVCB/mixture with no
                # single CID and no ladder entry. Report LOUDLY and drop it from
                # this uncertainty set (the resolvable members still run).
                logger.error("UNRESOLVABLE substance skipped (data-gap): %s", reason)
                warnings.warn(
                    f"Unresolvable substance accepted as data-gap and skipped: "
                    f"cas={cas!r} name={name!r}. Add a representative CAS to "
                    f"SURROGATE_CAS (survey/io.py) to include it.",
                    RuntimeWarning, stacklevel=2,
                )
                if data_gaps is not None:
                    data_gaps.append({"cas": cas, "name": name, "reason": reason})
                continue

            # Build SubstanceSpec
            # Use original CAS in the spec (for traceability), but note surrogate in name
            resolved_name = name or m.compound
            if used_surrogate:
                resolved_name = f"{resolved_name} [surrogate:{used_surrogate}]"

            # Get weight and exchangeable from spec
            # BACKWARD COMPATIBILITY: Default weight=1.0 if not specified
            weight = float(spec.get('weight', 1.0))

            # BACKWARD COMPATIBILITY: Default exchangeable=True if not specified
            # This preserves behavior for older scenarios without the 'exchangeable' field
            # Old behavior: all substances treated as exchangeable (mixture model)
            exchangeable_val = spec.get('exchangeable', True)
            if isinstance(exchangeable_val, bool):
                exchangeable = exchangeable_val
            else:
                # Handle 0/1, "true"/"false", etc.
                exchangeable = bool(int(exchangeable_val))  # 0/1 -> False/True

            sub_spec = SubstanceSpec.from_name(
                name=resolved_name,
                mass_g_mol=M,
                cas=cas or (m.CAS[0] if m.CAS else None),
            )
            # Set weight and exchangeable (from_name doesn't set these)
            sub_spec.weight = weight
            sub_spec.exchangeable = exchangeable
            substances.append(sub_spec)

        # A family whose EVERY member was unresolvable is a true failure, not a
        # silent empty result — surface it even under the accept-failure policy.
        if not substances and family_cfg.get('substances'):
            raise ValueError(
                "All substances in this family were unresolvable; "
                "no representative could be built. "
                f"data_gaps={data_gaps if data_gaps is not None else 'unrecorded'}"
            )
        return substances

    raise ValueError("Family config must have 'masses_g_mol' or 'substances'")


def load_scenario(path: Union[str, Path]) -> SurveyConfig:
    """
    Load scenario YAML and convert to SurveyConfig.

    Parameters
    ----------
    path : str or Path
        Path to scenario YAML file.

    Returns
    -------
    SurveyConfig
        Survey configuration.
    """
    data = load_yaml(path)

    physics = data.get('physics', {})
    priors = data.get('priors', {})
    solver = data.get('solver', {})
    fo_grid = solver.get('fo_grid', {})

    packaging = parse_packaging_spec(physics)

    # Patch P1 (2026-05-12): accept either `time_s` (1-step / monolayer)
    # or `step1_time_s` (two-step bilayer).  The latter is emitted by
    # convert_dbs_v3.py:298 when a reheating step is present, and was
    # silently dropped by the previous logic, causing every bilayer
    # 2-step manifest to model step 1 over the `PriorSpec` default
    # (max=10) instead of the YAML's true storage-time prior
    # (fixed in the v3 survey-engine consolidation).
    _step1 = priors.get('step1_time_s')
    _legacy = priors.get('time_s')
    if _step1 is not None and _legacy is not None and _step1 != _legacy:
        warnings.warn(
            f"Scenario {path} carries BOTH priors.step1_time_s and "
            "priors.time_s with different values; preferring step1_time_s. "
            "Reconcile the generator output.",
            RuntimeWarning, stacklevel=2,
        )
    time_prior_cfg = _step1 if _step1 is not None else (_legacy or {})
    time_prior = parse_prior_spec(time_prior_cfg, name='time')
    conc_prior = parse_prior_spec(priors.get('cp0_av', {}), name='concentration')

    return SurveyConfig(
        name=data.get('name', 'unnamed'),
        packaging=packaging,
        time_prior=time_prior,
        conc_prior=conc_prior,
        pdf_bins=int(solver.get('pdf_bins', 250)),
        n_fo=int(fo_grid.get('n_fo', 200)),
        fo_max_factor=float(fo_grid.get('fo_max_factor', 1.5)),
        fo_min_floor=float(fo_grid.get('fo_min_floor', 1e-15)),
        cache_dir=solver.get('cache_dir', '.survey_cache'),
        quantity=int(data.get('meta', {}).get('quantity', 1) or 1),
    )


def load_substances_from_scenario(
    path: Union[str, Path],
    *,
    skip_unresolvable: bool = True,
    data_gaps: Optional[List[Dict[str, Any]]] = None,
) -> List[SubstanceSpec]:
    """
    Load substances from scenario YAML.

    Parameters
    ----------
    path : str or Path
        Path to scenario YAML file.
    skip_unresolvable : bool, keyword-only, default True
        Forwarded to :func:`parse_substances` — accept unresolvable members
        as flagged data-gaps (True) or raise on the first one (False).
    data_gaps : list, keyword-only, optional
        Forwarded to :func:`parse_substances` to collect skipped members.

    Returns
    -------
    List[SubstanceSpec]
        List of substance specifications.
    """
    data = load_yaml(path)
    family = data.get('family', {})
    return parse_substances(
        family, skip_unresolvable=skip_unresolvable, data_gaps=data_gaps
    )


def save_results(
    path: Union[str, Path],
    pdf_bin_centers: np.ndarray,
    pdf: np.ndarray,
    cdf: np.ndarray,
    **arrays: np.ndarray,
) -> None:
    """
    Save survey results to NPZ file.

    Parameters
    ----------
    path : str or Path
        Output path.
    pdf_bin_centers : np.ndarray
        PDF bin centers.
    pdf : np.ndarray
        PDF values.
    cdf : np.ndarray
        CDF values.
    **arrays
        Additional arrays to save.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        pdf_bin_centers=pdf_bin_centers,
        pdf=pdf,
        cdf=cdf,
        **arrays,
    )


def save_manifest(
    path: Union[str, Path],
    config: SurveyConfig,
    substances: List[SubstanceSpec],
    fingerprint: str,
    **metadata: Any,
) -> None:
    """
    Save survey manifest (metadata) to JSON.

    Parameters
    ----------
    path : str or Path
        Output path.
    config : SurveyConfig
        Survey configuration.
    substances : List[SubstanceSpec]
        List of substances.
    fingerprint : str
        Survey fingerprint.
    **metadata
        Additional metadata.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    manifest = {
        'name': config.name,
        'fingerprint': fingerprint,
        'config': config.to_dict(),
        'substances': [s.to_dict() for s in substances],
        **metadata,
    }

    path.write_text(json.dumps(manifest, indent=2), encoding='utf-8')


# =============================================================================
# Multi-Component I/O (Phase 4)
# =============================================================================

def save_multicomponent_results(
    path: Union[str, Path],
    job_name: str,
    food_code: str,
    combined_cf: np.ndarray,
    combined_weights: np.ndarray,
    pdf_bin_centers: np.ndarray,
    pdf: np.ndarray,
    cdf: np.ndarray,
    component_results: Dict[str, Dict[str, Any]],
    statistics: Dict[str, float],
    **metadata: Any,
) -> None:
    """
    Save multi-component survey results to NPZ + JSON.

    Creates two files:
    - {path}.npz: NumPy arrays (CF samples, weights, PDF, CDF)
    - {path}_manifest.json: Metadata and per-component summary

    Parameters
    ----------
    path : str or Path
        Base output path (without extension).
    job_name : str
        Multi-component job name.
    food_code : str
        Food code identifier.
    combined_cf : np.ndarray
        Combined CF samples from tensor product.
    combined_weights : np.ndarray
        Combined weights from tensor product.
    pdf_bin_centers : np.ndarray
        PDF bin centers.
    pdf : np.ndarray
        Combined PDF values.
    cdf : np.ndarray
        Combined CDF values.
    component_results : Dict[str, Dict]
        Per-component results keyed by pack_code.
    statistics : Dict[str, float]
        Summary statistics (mean, std, quantiles).
    **metadata
        Additional metadata.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Save numerical arrays
    npz_path = path.with_suffix('.npz')
    np.savez_compressed(
        npz_path,
        combined_cf=combined_cf,
        combined_weights=combined_weights,
        pdf_bin_centers=pdf_bin_centers,
        pdf=pdf,
        cdf=cdf,
    )

    # Build manifest with component summaries
    component_summary = {}
    for pack_code, results in component_results.items():
        component_summary[pack_code] = {
            'polymer': results.get('polymer'),
            'n_samples': results.get('n_samples', len(results.get('CF_samples', []))),
            'substance_ids': results.get('substance_ids', []),
            'quantiles': results.get('quantiles', {}),
        }

    manifest = {
        'job_name': job_name,
        'food_code': food_code,
        'n_components': len(component_results),
        'pack_codes': list(component_results.keys()),
        'combined_n_samples': len(combined_cf),
        'statistics': statistics,
        'components': component_summary,
        **metadata,
    }

    json_path = path.parent / f"{path.stem}_manifest.json"
    json_path.write_text(json.dumps(manifest, indent=2), encoding='utf-8')


def load_multicomponent_results(path: Union[str, Path]) -> Dict[str, Any]:
    """
    Load multi-component survey results from NPZ + JSON.

    Parameters
    ----------
    path : str or Path
        Base path (with or without .npz extension).

    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - 'combined_cf': np.ndarray
        - 'combined_weights': np.ndarray
        - 'pdf_bin_centers': np.ndarray
        - 'pdf': np.ndarray
        - 'cdf': np.ndarray
        - 'manifest': dict with metadata
    """
    path = Path(path)
    if path.suffix != '.npz':
        npz_path = path.with_suffix('.npz')
    else:
        npz_path = path

    if not npz_path.exists():
        raise FileNotFoundError(f"Results file not found: {npz_path}")

    # Load arrays
    with np.load(npz_path) as data:
        results = {
            'combined_cf': data['combined_cf'],
            'combined_weights': data['combined_weights'],
            'pdf_bin_centers': data['pdf_bin_centers'],
            'pdf': data['pdf'],
            'cdf': data['cdf'],
        }

    # Load manifest
    json_path = npz_path.parent / f"{npz_path.stem}_manifest.json"
    if json_path.exists():
        results['manifest'] = json.loads(json_path.read_text(encoding='utf-8'))
    else:
        results['manifest'] = {}

    return results


def parse_multicomponent_scenario(data: Dict[str, Any]) -> MultiComponentJob:
    """
    Parse multi-component scenario from YAML dict.

    Expected format:
    ```yaml
    name: "E01_combined"
    food_code: "E01"
    food_volume_m3: 0.001
    food_simulant: "water"
    contact_temperature_degC: 25.0
    time_prior:
      triangular: {mode: 2592000, max: 7776000}
      grid: {nlow: 15, nhigh: 15}
    conc_prior:
      triangular: {mode: 50, max: 200}
      grid: {nlow: 15, nhigh: 15}
    components:
      - pack_code: "S1"
        polymer: "PET"
        thickness_m: 0.0001
        surface_area_m2: 0.06
        temperature_degC: 25.0
        substances:
          - {cas: "123-45-6", name: "SubstanceA", mass_g_mol: 150}
      - pack_code: "S2"
        polymer: "HDPE"
        thickness_m: 0.00015
        surface_area_m2: 0.04
        temperature_degC: 25.0
        substances:
          - {cas: "789-01-2", name: "SubstanceB", mass_g_mol: 200}
    ```

    Parameters
    ----------
    data : Dict[str, Any]
        Parsed YAML data.

    Returns
    -------
    MultiComponentJob
        Multi-component job specification.
    """
    # Parse priors
    time_prior_cfg = data.get('time_prior', {})
    conc_prior_cfg = data.get('conc_prior', {})

    time_prior = parse_prior_spec(time_prior_cfg, name='time')
    conc_prior = parse_prior_spec(conc_prior_cfg, name='concentration')

    # Parse components
    components = []
    for comp_cfg in data.get('components', []):
        # Parse substances for this component
        substances = []
        for sub_cfg in comp_cfg.get('substances', []):
            substances.append(SubstanceSpec(
                id=sub_cfg.get('cas', f"sub_{len(substances)}"),
                name=sub_cfg.get('name', ''),
                mass_g_mol=float(sub_cfg.get('mass_g_mol', 100)),
                weight=float(sub_cfg.get('weight', 1.0)),
            ))

        components.append(ComponentSpec(
            pack_code=comp_cfg['pack_code'],
            polymer=comp_cfg['polymer'],
            thickness_m=float(comp_cfg.get('thickness_m', 100e-6)),
            surface_area_m2=float(comp_cfg.get('surface_area_m2', 0.06)),
            temperature_degC=float(comp_cfg.get('temperature_degC', 25.0)),
            substances=substances,
        ))

    return MultiComponentJob(
        name=data.get('name', 'unnamed'),
        food_code=data.get('food_code', ''),
        components=components,
        food_volume_m3=float(data.get('food_volume_m3', 0.001)),
        food_simulant=data.get('food_simulant', 'water'),
        contact_temperature_degC=float(data.get('contact_temperature_degC', 25.0)),
        time_prior=time_prior,
        conc_prior=conc_prior,
        h_m_s=float(data.get('h_m_s', 1e-7)),
        cf0=float(data.get('cf0', 0.0)),
    )


def load_multicomponent_scenario(path: Union[str, Path]) -> MultiComponentJob:
    """
    Load multi-component scenario from YAML file.

    Parameters
    ----------
    path : str or Path
        Path to YAML file.

    Returns
    -------
    MultiComponentJob
        Multi-component job specification.
    """
    data = load_yaml(path)
    return parse_multicomponent_scenario(data)


def save_multicomponent_scenario(
    path: Union[str, Path],
    job: MultiComponentJob,
) -> None:
    """
    Save multi-component job specification to YAML.

    Parameters
    ----------
    path : str or Path
        Output path.
    job : MultiComponentJob
        Job specification.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    data = {
        'name': job.name,
        'food_code': job.food_code,
        'food_volume_m3': job.food_volume_m3,
        'food_simulant': job.food_simulant,
        'contact_temperature_degC': job.contact_temperature_degC,
        'h_m_s': job.h_m_s,
        'cf0': job.cf0,
        'time_prior': {
            'triangular': {
                'mode': job.time_prior.mode,
                'max': job.time_prior.max_val,
            },
            'grid': {
                'nlow': job.time_prior.n_low,
                'nhigh': job.time_prior.n_high,
            },
        },
        'conc_prior': {
            'triangular': {
                'mode': job.conc_prior.mode,
                'max': job.conc_prior.max_val,
            },
            'grid': {
                'nlow': job.conc_prior.n_low,
                'nhigh': job.conc_prior.n_high,
            },
        },
        'components': [
            {
                'pack_code': comp.pack_code,
                'polymer': comp.polymer,
                'thickness_m': comp.thickness_m,
                'surface_area_m2': comp.surface_area_m2,
                'temperature_degC': comp.temperature_degC,
                'substances': [
                    {
                        'cas': s.id,
                        'name': s.name,
                        'mass_g_mol': s.mass_g_mol,
                        'weight': s.weight,
                    }
                    for s in comp.substances
                ],
            }
            for comp in job.components
        ],
    }

    with open(path, 'w', encoding='utf-8') as f:
        yaml.dump(data, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
