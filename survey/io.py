"""
survey/io.py — Input/Output Functions for Survey
================================================

Handles scenario loading from YAML and result dumping to NPZ/JSON.

@project: SFPPy/INSERM — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

import json
from pathlib import Path
from typing import Dict, Any, Union, List, Optional

import numpy as np
import yaml

from survey.models import (
    LayerSpec,
    PackagingSpec,
    SubstanceSpec,
    PriorSpec,
    SurveyConfig,
)


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


def parse_layer_specs(physics: Dict[str, Any]) -> List[LayerSpec]:
    """
    Parse layer specifications from physics config.

    Supports both monolayer and multilayer formats.

    Parameters
    ----------
    physics : dict
        Physics section from scenario YAML.

    Returns
    -------
    List[LayerSpec]
        List of layer specifications.
    """
    # Check for monolayer format
    if 'monolayer' in physics:
        mono = physics['monolayer']
        return [LayerSpec(
            polymer=mono.get('polymer', 'LDPE'),
            thickness_m=float(mono.get('thickness_m', 100e-6)),
            C0=float(mono.get('C0', 0.0)),
            temperature_degC=float(mono.get('temperature_degC', 40.0)),
        )]

    # Check for multilayer format
    if 'multilayer' in physics:
        multi = physics['multilayer']
        layers_cfg = multi.get('layers', [])
        T = float(multi.get('temperature_degC', 40.0))
        return [
            LayerSpec(
                polymer=layer.get('polymer', 'LDPE'),
                thickness_m=float(layer.get('thickness_m', 100e-6)),
                C0=float(layer.get('C0', 0.0)),
                temperature_degC=T,
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
    layers = parse_layer_specs(physics)
    iface = physics.get('interface', {})

    return PackagingSpec(
        layers=layers,
        h_m_s=float(iface.get('h_m_s', 1e-7)),
        surface_area_m2=float(iface.get('surface_area_m2', 0.06)),
        food_volume_m3=float(iface.get('food_volume_m3', 0.001)),
        contact_temperature_degC=float(iface.get('contact_temperature_degC', 40.0)),
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
        mode=float(tri.get('mode', 1.0)),
        max_val=float(tri.get('max', 10.0)),
        n_low=int(grid.get('nlow', 15)),
        n_high=int(grid.get('nhigh', 15)),
        name=name,
    )


def parse_substances(family_cfg: Dict[str, Any]) -> List[SubstanceSpec]:
    """
    Parse substance specifications from family config.

    Supports:
    1. Direct masses: family.masses_g_mol = [80, 120, ...]
    2. Named substances: family.substances = [{name: "...", cas: "..."}, ...]

    Parameters
    ----------
    family_cfg : dict
        Family section from scenario YAML.

    Returns
    -------
    List[SubstanceSpec]
        List of substance specifications.
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
        from patankar.loadpubchem import migrant

        substances = []
        for i, spec in enumerate(family_cfg['substances']):
            name = spec.get('name')
            cas = spec.get('cas')

            if name is None and cas is None:
                raise ValueError(f"Substance spec must have 'name' or 'cas': {spec}")

            # Fetch molar mass from PubChem
            # Try CAS first (more reliable), then fall back to name
            try:
                m = migrant(name=cas or name)
                M = float(m.M)
            except Exception as e:
                # If CAS failed and we have a name, try name as fallback
                if cas and name:
                    try:
                        m = migrant(name=name)
                        M = float(m.M)
                    except Exception:
                        raise ValueError(f"Failed to resolve substance '{cas}' (CAS) or '{name}' (name): {e}")
                else:
                    raise ValueError(f"Failed to resolve substance '{cas or name}': {e}")

            substances.append(SubstanceSpec.from_name(
                name=name or m.compound,
                mass_g_mol=M,
                cas=cas or (m.CAS[0] if m.CAS else None),
            ))

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
    time_prior = parse_prior_spec(priors.get('time_s', {}), name='time')
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
    )


def load_substances_from_scenario(path: Union[str, Path]) -> List[SubstanceSpec]:
    """
    Load substances from scenario YAML.

    Parameters
    ----------
    path : str or Path
        Path to scenario YAML file.

    Returns
    -------
    List[SubstanceSpec]
        List of substance specifications.
    """
    data = load_yaml(path)
    family = data.get('family', {})
    return parse_substances(family)


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
