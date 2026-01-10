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
    ComponentSpec,
    MultiComponentJob,
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
