"""
survey/fingerprints.py — Canonical Hashing and Fingerprints
==========================================================

Provides deterministic fingerprinting for survey components:
- Physics fingerprint: (packaging, substance, solver params)
- Probability fingerprint: (priors, substance set)

@project: SFPPy/INSERM — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

import hashlib
import json
from typing import Any, Dict, List
from dataclasses import asdict

from survey.models import (
    LayerSpec,
    PackagingSpec,
    SubstanceSpec,
    PriorSpec,
    SurveyConfig,
)


def canonical_json(obj: Any, float_precision: int = 10) -> str:
    """
    Convert object to canonical JSON string.

    Ensures deterministic serialization:
    - Sorted keys
    - Compact separators
    - Stable float formatting

    Parameters
    ----------
    obj : Any
        Object to serialize.
    float_precision : int
        Number of decimal places for floats.

    Returns
    -------
    str
        Canonical JSON string.
    """
    def normalize(v):
        if isinstance(v, float):
            return round(v, float_precision)
        if isinstance(v, dict):
            return {k: normalize(val) for k, val in sorted(v.items())}
        if isinstance(v, (list, tuple)):
            return [normalize(item) for item in v]
        return v

    normalized = normalize(obj)
    return json.dumps(normalized, sort_keys=True, separators=(",", ":"))


def stable_hash(data: str) -> str:
    """
    Compute SHA-256 hash of string data.

    Parameters
    ----------
    data : str
        String to hash.

    Returns
    -------
    str
        Hexadecimal hash string.
    """
    return hashlib.sha256(data.encode("utf-8")).hexdigest()


def fingerprint_physics(
    polymer: str,
    mass_g_mol: float,
    D: float,
    k: float,
    k0: float,
    lP_m: float,
    Fo_max: float,
    n_fo: int,
    h: float,
    surface_area: float,
    food_volume: float,
    contact_temperature_degC: float,
    CF0: float,
) -> str:
    """
    Compute physics fingerprint for master curve caching.

    This fingerprint uniquely identifies a master curve computation.
    Changes to any parameter require recomputation.

    Parameters
    ----------
    polymer : str
        Polymer identifier.
    mass_g_mol : float
        Substance molar mass.
    D : float
        Diffusivity (m²/s).
    k : float
        Partition coefficient.
    k0 : float
        k0 parameter.
    lP_m : float
        Layer thickness (m).
    Fo_max : float
        Maximum Fourier number.
    n_fo : int
        Number of Fo grid points.
    h : float
        Mass transfer coefficient (m/s).
    surface_area : float
        Contact surface area (m²).
    food_volume : float
        Food volume (m³).
    contact_temperature_degC : float
        Contact temperature (°C).
    CF0 : float
        Initial food concentration.

    Returns
    -------
    str
        64-character hexadecimal fingerprint.
    """
    payload = {
        "polymer": polymer,
        "mass_g_mol": mass_g_mol,
        "D": D,
        "k": k,
        "k0": k0,
        "lP_m": lP_m,
        "Fo_max": Fo_max,
        "n_fo": n_fo,
        "h": h,
        "surface_area": surface_area,
        "food_volume": food_volume,
        "contact_temperature_degC": contact_temperature_degC,
        "CF0": CF0,
    }
    return stable_hash(canonical_json(payload))


def fingerprint_substance_set(substances: List[SubstanceSpec]) -> str:
    """
    Compute fingerprint for a set of substances.

    Parameters
    ----------
    substances : List[SubstanceSpec]
        List of substances.

    Returns
    -------
    str
        Fingerprint of the substance set.
    """
    # Sort by canonical ID for determinism
    sorted_subs = sorted(substances, key=lambda s: s.canonical_id())
    payload = [
        {
            "id": s.id,
            "mass_g_mol": s.mass_g_mol,
            "D": s.D,
            "k": s.k,
            "k0": s.k0,
        }
        for s in sorted_subs
    ]
    return stable_hash(canonical_json(payload))


def fingerprint_prior(prior: PriorSpec) -> str:
    """
    Compute fingerprint for a prior specification.

    Parameters
    ----------
    prior : PriorSpec
        Prior specification.

    Returns
    -------
    str
        Fingerprint.
    """
    return stable_hash(canonical_json(prior.to_dict()))


def fingerprint_packaging(packaging: PackagingSpec) -> str:
    """
    Compute fingerprint for packaging specification.

    Parameters
    ----------
    packaging : PackagingSpec
        Packaging specification.

    Returns
    -------
    str
        Fingerprint.
    """
    return stable_hash(canonical_json(packaging.to_dict()))


def fingerprint_probability(
    substances: List[SubstanceSpec],
    packaging: PackagingSpec,
    time_prior: PriorSpec,
    conc_prior: PriorSpec,
    i_ref: int = 0,
) -> str:
    """
    Compute probability fingerprint for PDF caching.

    This fingerprint identifies when a PDF needs recomputation.
    Changes to substances, priors, or reference layer require
    PDF recomputation.

    Parameters
    ----------
    substances : List[SubstanceSpec]
        List of substances with inferred parameters.
    packaging : PackagingSpec
        Packaging specification.
    time_prior : PriorSpec
        Time prior.
    conc_prior : PriorSpec
        Concentration prior.
    i_ref : int
        Reference layer index (for multilayer).

    Returns
    -------
    str
        Fingerprint for probability cache.
    """
    payload = {
        "substances": fingerprint_substance_set(substances),
        "packaging": fingerprint_packaging(packaging),
        "time_prior": fingerprint_prior(time_prior),
        "conc_prior": fingerprint_prior(conc_prior),
        "i_ref": i_ref,
    }
    return stable_hash(canonical_json(payload))


def fingerprint_survey_config(config: SurveyConfig) -> str:
    """
    Compute full fingerprint for survey configuration.

    Parameters
    ----------
    config : SurveyConfig
        Full survey configuration.

    Returns
    -------
    str
        Survey fingerprint.
    """
    return stable_hash(canonical_json(config.to_dict()))
