"""
survey/aggregation.py — Multi-Component Tensor Product Aggregation
====================================================================

Combines migration contributions from multiple packaging components
using deterministic tensor product expansion (no Monte Carlo, no FFT).

Theory
------
For independent random variables X₁, X₂ with discrete probability masses:
    X_total = X₁ + X₂

The PDF of X_total is computed via tensor product expansion:
    CF_total[i,j] = CF₁[i] + CF₂[j]
    W_total[i,j] = W₁[i] × W₂[j]

This is the deterministic equivalent of convolution, consistent with
survey's finite-difference quadrature approach.

@project: SFPPy/INSERM — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

import numpy as np
from typing import Dict, Tuple, List, Optional, Any


def combine_tensors(
    cf1: np.ndarray,
    w1: np.ndarray,
    cf2: np.ndarray,
    w2: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Combine two CF tensors via tensor product (sum of independent variables).

    Uses outer product expansion, NOT FFT convolution.
    Consistent with survey's deterministic discretization approach.

    Mathematical basis:
        CF_total[i,j] = CF₁[i] + CF₂[j]   (addition of values)
        W_total[i,j] = W₁[i] × W₂[j]      (product of probabilities)

    Parameters
    ----------
    cf1 : np.ndarray
        Flattened CF samples from component 1 (1D array).
    w1 : np.ndarray
        Corresponding weights for cf1 (sum to 1).
    cf2 : np.ndarray
        Flattened CF samples from component 2 (1D array).
    w2 : np.ndarray
        Corresponding weights for cf2 (sum to 1).

    Returns
    -------
    cf_combined : np.ndarray
        Combined CF samples (shape: len(cf1) × len(cf2), flattened).
    w_combined : np.ndarray
        Combined weights (outer product, renormalized to sum to 1).

    Examples
    --------
    >>> cf1 = np.array([0.1, 0.2, 0.3])
    >>> w1 = np.array([0.2, 0.5, 0.3])
    >>> cf2 = np.array([0.05, 0.15])
    >>> w2 = np.array([0.4, 0.6])
    >>> cf_combined, w_combined = combine_tensors(cf1, w1, cf2, w2)
    >>> cf_combined.shape
    (6,)
    >>> np.isclose(w_combined.sum(), 1.0)
    True
    """
    # Ensure 1D arrays
    cf1 = np.asarray(cf1).ravel()
    cf2 = np.asarray(cf2).ravel()
    w1 = np.asarray(w1, dtype=float).ravel()
    w2 = np.asarray(w2, dtype=float).ravel()

    # Tensor product: CF_total[i,j] = CF₁[i] + CF₂[j]
    cf_combined = np.add.outer(cf1, cf2).ravel()

    # Weight product: W_total[i,j] = W₁[i] × W₂[j]
    w_combined = np.outer(w1, w2).ravel()

    # Renormalize (should already be ~1 if inputs are normalized)
    w_sum = w_combined.sum()
    if w_sum > 0:
        w_combined /= w_sum

    return cf_combined, w_combined


def combine_multiple_tensors(
    tensors: List[Tuple[np.ndarray, np.ndarray]]
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Combine multiple CF tensors via sequential tensor product.

    For n components: CF_total = CF₁ + CF₂ + ... + CFₙ
    Computed by chaining pairwise combinations.

    Parameters
    ----------
    tensors : List[Tuple[np.ndarray, np.ndarray]]
        List of (CF_samples, weights) tuples from each component.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Combined (CF_samples, weights).

    Raises
    ------
    ValueError
        If tensors list is empty.

    Notes
    -----
    Computational complexity:
    - 2 components: O(n₁ × n₂)
    - 3 components: O(n₁ × n₂ × n₃)
    - etc.

    For typical 15×15 = 225 samples per component:
    - 2 components: 50,625 samples
    - 3 components: ~11M samples
    """
    if not tensors:
        raise ValueError("Cannot combine empty list of tensors")

    if len(tensors) == 1:
        return tensors[0]

    cf, w = tensors[0]
    for cf2, w2 in tensors[1:]:
        cf, w = combine_tensors(cf, w, cf2, w2)

    return cf, w


def aggregate_components(
    component_results: Dict[str, Dict[str, Any]],
    substance_key: str = 'substance_ids'
) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Aggregate migration from multiple components for each substance.

    For each substance present in ANY component:
    - If in one component: use that tensor directly
    - If in multiple: combine via tensor product expansion

    Parameters
    ----------
    component_results : Dict[str, dict]
        Results from each component keyed by pack_code.
        Each dict must contain:
        - 'CF_samples': np.ndarray of CF values
        - 'weights': np.ndarray of corresponding weights
        - substance_key: list of substance IDs

    substance_key : str
        Key in results dict containing substance IDs list.
        Default: 'substance_ids'

    Returns
    -------
    Dict[str, Dict[str, np.ndarray]]
        Aggregated results keyed by substance_id.
        Each dict contains:
        - 'CF_samples': combined CF samples
        - 'weights': combined weights
        - 'source_components': list of pack_codes contributing

    Examples
    --------
    >>> # Two components, substance A in both, substance B only in S1
    >>> results = {
    ...     'S1': {'CF_samples': np.array([0.1, 0.2]), 'weights': np.array([0.5, 0.5]),
    ...            'substance_ids': ['A', 'B']},
    ...     'S2': {'CF_samples': np.array([0.05, 0.15]), 'weights': np.array([0.6, 0.4]),
    ...            'substance_ids': ['A']}
    ... }
    >>> agg = aggregate_components(results)
    >>> 'A' in agg and 'B' in agg
    True
    >>> len(agg['A']['source_components'])  # A from both S1 and S2
    2
    >>> len(agg['B']['source_components'])  # B only from S1
    1
    """
    # Collect all unique substances across components
    all_substances = set()
    for results in component_results.values():
        substance_ids = results.get(substance_key, [])
        if substance_ids:
            all_substances.update(substance_ids)

    aggregated = {}
    for substance_id in all_substances:
        # Collect (CF, weights) from each component where substance is present
        tensors = []
        source_components = []

        for pack_code, results in component_results.items():
            substance_ids = results.get(substance_key, [])
            if substance_id in substance_ids:
                tensors.append((
                    np.asarray(results['CF_samples']),
                    np.asarray(results['weights'])
                ))
                source_components.append(pack_code)

        # Combine via tensor product
        if tensors:
            cf, w = combine_multiple_tensors(tensors)
            aggregated[substance_id] = {
                'CF_samples': cf,
                'weights': w,
                'source_components': source_components
            }

    return aggregated


def aggregate_family_weighted(
    substance_results: Dict[str, Dict[str, np.ndarray]],
    substance_weights: Dict[str, float]
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Aggregate multiple substances into a family PDF using occurrence weights.

    Each substance contributes to the family PDF proportionally to its
    occurrence weight (e.g., 0.5, 1.0, 2.0 from database).

    Mathematical basis:
    For substances s₁, s₂, ... with weights ω₁, ω₂, ...:
        W_family[k] = Σᵢ (ωᵢ / Σω) × W_sᵢ[k]

    Parameters
    ----------
    substance_results : Dict[str, Dict[str, np.ndarray]]
        Results keyed by substance_id, each containing
        'CF_samples' and 'weights'.

    substance_weights : Dict[str, float]
        Occurrence weights for each substance.
        Missing substances are assigned weight 1.0.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (CF_samples_combined, weights_combined) for the family.

    Notes
    -----
    This is different from tensor product aggregation.
    Here we're computing a weighted mixture of PDFs, not a sum of RVs.

    For family exposure, substances are alternatives (OR), not additive (AND).
    The family PDF represents "exposure given any substance in family".
    """
    if not substance_results:
        raise ValueError("Cannot aggregate empty substance results")

    # Normalize occurrence weights
    total_weight = sum(
        substance_weights.get(sid, 1.0)
        for sid in substance_results.keys()
    )
    if total_weight <= 0:
        total_weight = len(substance_results)  # Fallback to equiprobable

    # Find maximum sample length (for padding if needed)
    max_len = max(len(r['CF_samples']) for r in substance_results.values())

    # Accumulate weighted contributions
    cf_all = []
    w_all = []

    for substance_id, results in substance_results.items():
        cf = np.asarray(results['CF_samples'])
        w = np.asarray(results['weights'])

        # Apply occurrence weight
        occurrence_weight = substance_weights.get(substance_id, 1.0) / total_weight
        w_weighted = w * occurrence_weight

        cf_all.append(cf)
        w_all.append(w_weighted)

    # Concatenate all samples (mixture model)
    cf_combined = np.concatenate(cf_all)
    w_combined = np.concatenate(w_all)

    # Renormalize
    w_sum = w_combined.sum()
    if w_sum > 0:
        w_combined /= w_sum

    return cf_combined, w_combined


def compute_pdf_from_samples(
    cf_samples: np.ndarray,
    weights: np.ndarray,
    n_bins: int = 250,
    cf_min: Optional[float] = None,
    cf_max: Optional[float] = None
) -> Dict[str, np.ndarray]:
    """
    Compute PDF and CDF from weighted samples.

    Parameters
    ----------
    cf_samples : np.ndarray
        CF values (1D array).
    weights : np.ndarray
        Corresponding weights (sum to 1).
    n_bins : int
        Number of histogram bins.
    cf_min : float, optional
        Minimum bin edge. Defaults to 0.
    cf_max : float, optional
        Maximum bin edge. Defaults to max(cf_samples) * 1.1.

    Returns
    -------
    Dict[str, np.ndarray]
        Dictionary containing:
        - 'pdf': probability density values
        - 'cdf': cumulative distribution values
        - 'bin_centers': bin center values
        - 'bin_edges': bin edge values
    """
    cf_samples = np.asarray(cf_samples)
    weights = np.asarray(weights, dtype=float)

    # Handle edge case: all zeros
    if cf_samples.max() == 0:
        bin_edges = np.linspace(0, 1, n_bins + 1)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        pdf = np.zeros(n_bins)
        pdf[0] = 1.0 / (bin_edges[1] - bin_edges[0])
        cdf = np.ones(n_bins)
        return {
            'pdf': pdf,
            'cdf': cdf,
            'bin_centers': bin_centers,
            'bin_edges': bin_edges
        }

    # Set bin range
    if cf_min is None:
        cf_min = 0.0
    if cf_max is None:
        cf_max = cf_samples.max() * 1.1

    bin_edges = np.linspace(cf_min, cf_max, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width = bin_edges[1] - bin_edges[0]

    # Weighted histogram
    pdf_counts, _ = np.histogram(cf_samples, bins=bin_edges, weights=weights)

    # Normalize to PDF (integral = 1)
    pdf = pdf_counts / (bin_width * pdf_counts.sum()) if pdf_counts.sum() > 0 else pdf_counts

    # Compute CDF
    cdf = np.cumsum(pdf) * bin_width

    return {
        'pdf': pdf,
        'cdf': cdf,
        'bin_centers': bin_centers,
        'bin_edges': bin_edges
    }


def quantile_from_samples(
    cf_samples: np.ndarray,
    weights: np.ndarray,
    q: float
) -> float:
    """
    Compute quantile from weighted samples (exact empirical CDF).

    This function computes quantiles directly from the weighted samples
    without binning, providing exact empirical CDF values.

    Parameters
    ----------
    cf_samples : np.ndarray
        CF values (1D array).
    weights : np.ndarray
        Corresponding weights (sum to 1).
    q : float
        Quantile (0 to 1).

    Returns
    -------
    float
        Quantile value.
    """
    cf_samples = np.asarray(cf_samples)
    weights = np.asarray(weights, dtype=float)

    # Sort samples and weights
    sorted_indices = np.argsort(cf_samples)
    sorted_cf = cf_samples[sorted_indices]
    sorted_w = weights[sorted_indices]

    # Compute cumulative weights
    cumsum = np.cumsum(sorted_w)
    cumsum /= cumsum[-1]  # Normalize

    # Find quantile
    idx = np.searchsorted(cumsum, q)
    if idx >= len(sorted_cf):
        return sorted_cf[-1]
    return sorted_cf[idx]


def quantile_from_cdf(
    bin_centers: np.ndarray,
    cdf: np.ndarray,
    q: float
) -> float:
    """
    Compute quantile from integrated CDF (PDF histogram method).

    This function extracts quantiles from a CDF that was computed
    by numerical integration of a PDF histogram. Use this for
    consistency with stored NPZ outputs.

    Parameters
    ----------
    bin_centers : np.ndarray
        Bin center values from PDF histogram (1D array).
    cdf : np.ndarray
        Cumulative distribution function values at bin centers.
    q : float
        Quantile (0 to 1).

    Returns
    -------
    float
        Quantile value interpolated from CDF.

    Notes
    -----
    This provides values consistent with the PDF/CDF stored in NPZ files.
    For risk assessment, use compute_risk_percentiles() which uses the
    exact empirical CDF from samples.
    """
    bin_centers = np.asarray(bin_centers)
    cdf = np.asarray(cdf)

    # Find index where CDF crosses the quantile
    idx = np.searchsorted(cdf, q)
    if idx >= len(bin_centers):
        return float(bin_centers[-1])
    if idx == 0:
        return float(bin_centers[0])

    # Linear interpolation for more accurate value
    cdf_lo, cdf_hi = cdf[idx - 1], cdf[idx]
    x_lo, x_hi = bin_centers[idx - 1], bin_centers[idx]

    if cdf_hi - cdf_lo > 1e-10:
        frac = (q - cdf_lo) / (cdf_hi - cdf_lo)
        return float(x_lo + frac * (x_hi - x_lo))
    return float(x_hi)


def compute_statistics(
    cf_samples: np.ndarray,
    weights: np.ndarray
) -> Dict[str, float]:
    """
    Compute weighted statistics from samples.

    Parameters
    ----------
    cf_samples : np.ndarray
        CF values (1D array).
    weights : np.ndarray
        Corresponding weights.

    Returns
    -------
    Dict[str, float]
        Statistics: mean, std, q50, q95, q99, max.
    """
    cf_samples = np.asarray(cf_samples)
    weights = np.asarray(weights, dtype=float)

    # Normalize weights
    w_sum = weights.sum()
    if w_sum > 0:
        w = weights / w_sum
    else:
        w = np.ones_like(weights) / len(weights)

    mean = float(np.average(cf_samples, weights=w))
    variance = float(np.average((cf_samples - mean)**2, weights=w))
    std = float(np.sqrt(variance))

    return {
        'mean': mean,
        'std': std,
        'q50': float(quantile_from_samples(cf_samples, w, 0.50)),
        'q95': float(quantile_from_samples(cf_samples, w, 0.95)),
        'q99': float(quantile_from_samples(cf_samples, w, 0.99)),
        'max': float(cf_samples.max())
    }


# Standard percentiles for risk assessment (EFSA, FDA, regulatory contexts)
RISK_PERCENTILES = [50, 75, 90, 95, 97.5, 99]


def compute_risk_percentiles(
    cf_samples: np.ndarray,
    weights: np.ndarray
) -> Dict[str, float]:
    """
    Compute risk-relevant percentiles for regulatory assessment.

    These are the standard percentiles used by risk scientists:
    - P50: Median (central tendency)
    - P75: Third quartile
    - P90: High exposure
    - P95: Reasonable worst case (commonly used threshold)
    - P97.5: Used in some EFSA/regulatory frameworks
    - P99: Extreme exposure

    Parameters
    ----------
    cf_samples : np.ndarray
        CF values (1D array).
    weights : np.ndarray
        Corresponding weights.

    Returns
    -------
    Dict[str, float]
        Risk percentiles with keys: 'p50', 'p75', 'p90', 'p95', 'p97_5', 'p99'

    Notes
    -----
    These percentiles are independent of the resolution parameter in
    compute_percentiles(). They provide the exact values needed for
    regulatory reporting.
    """
    cf_samples = np.asarray(cf_samples)
    weights = np.asarray(weights, dtype=float)

    # Normalize weights
    w_sum = weights.sum()
    if w_sum > 0:
        w = weights / w_sum
    else:
        w = np.ones_like(weights) / len(weights)

    return {
        'p50': float(quantile_from_samples(cf_samples, w, 0.50)),
        'p75': float(quantile_from_samples(cf_samples, w, 0.75)),
        'p90': float(quantile_from_samples(cf_samples, w, 0.90)),
        'p95': float(quantile_from_samples(cf_samples, w, 0.95)),
        'p97_5': float(quantile_from_samples(cf_samples, w, 0.975)),
        'p99': float(quantile_from_samples(cf_samples, w, 0.99)),
    }


def compute_percentiles(
    cf_samples: np.ndarray,
    weights: np.ndarray,
    resolution: int = 50
) -> Dict[str, Any]:
    """
    Compute CDF percentiles at configurable resolution.

    Parameters
    ----------
    cf_samples : np.ndarray
        CF values (1D array).
    weights : np.ndarray
        Corresponding weights.
    resolution : int
        Number of percentiles (default 50 → P2, P4, ..., P100).

    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - 'levels': np.ndarray of percentile levels (e.g., [2, 4, 6, ..., 100])
        - 'values': np.ndarray of corresponding CF values
        - 'percentiles': Dict[str, float] with keys 'p02', 'p04', etc.

    Examples
    --------
    >>> cf = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    >>> w = np.ones(5) / 5
    >>> result = compute_percentiles(cf, w, resolution=50)
    >>> len(result['values'])
    50
    >>> 'p50' in result['percentiles']
    True
    """
    cf_samples = np.asarray(cf_samples)
    weights = np.asarray(weights, dtype=float)

    # Normalize weights
    w_sum = weights.sum()
    if w_sum > 0:
        w = weights / w_sum
    else:
        w = np.ones_like(weights) / len(weights)

    # Generate percentile levels (e.g., 2, 4, 6, ..., 100 for resolution=50)
    step = 100 / resolution
    levels = np.arange(step, 100 + step/2, step)  # Include 100
    levels = np.minimum(levels, 100)  # Ensure max is 100

    # Compute quantiles for each level
    values = np.array([
        quantile_from_samples(cf_samples, w, level / 100.0)
        for level in levels
    ])

    # Build percentile dictionary with formatted keys
    percentiles = {}
    for level, value in zip(levels, values):
        key = f"p{int(level):02d}"
        percentiles[key] = float(value)

    return {
        'levels': levels,
        'values': values,
        'percentiles': percentiles
    }


def aggregate_packaging_components(
    component_results: List[Dict[str, Any]],
    resolution: int = 50
) -> Dict[str, Any]:
    """
    Aggregate CDFs from multiple packaging components via tensor product.

    This is the main entry point for multi-component food packaging aggregation.
    Combines migration from independent components (e.g., PET bottle + HDPE cap).

    Parameters
    ----------
    component_results : List[Dict[str, Any]]
        Results from each packaging component. Each dict must contain:
        - 'cf_samples': np.ndarray of CF values
        - 'weights': np.ndarray of corresponding weights
        - 'pack_code': str identifying the component
        Optional:
        - 'polymer': str polymer type
        - 'family': str tech function family

    resolution : int
        Number of percentiles for output CDF (default 50).

    Returns
    -------
    Dict[str, Any]
        Aggregated results containing:
        - 'cf_samples': combined CF samples (tensor product)
        - 'weights': combined weights
        - 'pdf': probability density function
        - 'cdf': cumulative distribution function
        - 'bin_centers': PDF/CDF x-axis values
        - 'percentiles': Dict with 'levels', 'values', 'percentiles'
        - 'statistics': mean, std, q50, q95, q99, max
        - 'components': list of input component metadata

    Examples
    --------
    >>> # Two-component packaging: PET bottle body + HDPE cap
    >>> comp1 = {'cf_samples': np.array([0.1, 0.2]), 'weights': np.array([0.5, 0.5]),
    ...          'pack_code': 'A12_S1', 'polymer': 'PET'}
    >>> comp2 = {'cf_samples': np.array([0.05, 0.15]), 'weights': np.array([0.6, 0.4]),
    ...          'pack_code': 'A12_S2', 'polymer': 'HDPE'}
    >>> result = aggregate_packaging_components([comp1, comp2], resolution=50)
    >>> result['cf_samples'].shape[0]  # 2 × 2 = 4 combined samples
    4
    >>> len(result['percentiles']['values'])
    50
    """
    if not component_results:
        raise ValueError("Cannot aggregate empty component list")

    # Extract tensors
    tensors = []
    components_meta = []
    for comp in component_results:
        cf = np.asarray(comp['cf_samples'])
        w = np.asarray(comp['weights'])
        tensors.append((cf, w))
        components_meta.append({
            'pack_code': comp.get('pack_code', 'unknown'),
            'polymer': comp.get('polymer', 'unknown'),
            'family': comp.get('family', 'unknown'),
            'n_samples': len(cf)
        })

    # Combine via tensor product
    cf_combined, w_combined = combine_multiple_tensors(tensors)

    # Compute PDF/CDF
    pdf_result = compute_pdf_from_samples(cf_combined, w_combined)

    # Compute percentiles at requested resolution
    percentile_result = compute_percentiles(cf_combined, w_combined, resolution)

    # Compute statistics
    stats = compute_statistics(cf_combined, w_combined)

    return {
        'cf_samples': cf_combined,
        'weights': w_combined,
        'pdf': pdf_result['pdf'],
        'cdf': pdf_result['cdf'],
        'bin_centers': pdf_result['bin_centers'],
        'bin_edges': pdf_result['bin_edges'],
        'percentiles': percentile_result,
        'statistics': stats,
        'components': components_meta,
        'n_components': len(component_results),
        'n_combined_samples': len(cf_combined)
    }


def aggregate_food_exposure(
    packaging_families: Dict[str, List[Dict[str, Any]]],
    resolution: int = 50
) -> Dict[str, Dict[str, Any]]:
    """
    Aggregate exposure for a food item across all tech function families.

    For each family (antioxidant, plasticizer, etc.), aggregates contributions
    from all packaging components, then returns family-level CDFs.

    Parameters
    ----------
    packaging_families : Dict[str, List[Dict[str, Any]]]
        Results organized by family name. Each family contains a list of
        component results (same structure as aggregate_packaging_components input).

    resolution : int
        Number of percentiles for output CDFs (default 50).

    Returns
    -------
    Dict[str, Dict[str, Any]]
        Aggregated results keyed by family name.
        Each family dict contains the same structure as aggregate_packaging_components output.

    Examples
    --------
    >>> # Food A12 with antioxidants from PET and HDPE components
    >>> families = {
    ...     'antioxidant': [
    ...         {'cf_samples': np.array([0.1, 0.2]), 'weights': np.array([0.5, 0.5]),
    ...          'pack_code': 'A12_S1', 'polymer': 'PET'},
    ...         {'cf_samples': np.array([0.3, 0.4]), 'weights': np.array([0.4, 0.6]),
    ...          'pack_code': 'A12_S2', 'polymer': 'HDPE'}
    ...     ],
    ...     'plasticizer': [
    ...         {'cf_samples': np.array([0.05, 0.1]), 'weights': np.array([0.7, 0.3]),
    ...          'pack_code': 'A12_S1', 'polymer': 'PET'}
    ...     ]
    ... }
    >>> result = aggregate_food_exposure(families, resolution=50)
    >>> 'antioxidant' in result and 'plasticizer' in result
    True
    >>> result['antioxidant']['n_components']
    2
    """
    aggregated = {}

    for family_name, component_list in packaging_families.items():
        if not component_list:
            continue

        if len(component_list) == 1:
            # Single component: just compute percentiles, no tensor product needed
            comp = component_list[0]
            cf = np.asarray(comp['cf_samples'])
            w = np.asarray(comp['weights'])

            pdf_result = compute_pdf_from_samples(cf, w)
            percentile_result = compute_percentiles(cf, w, resolution)
            stats = compute_statistics(cf, w)

            aggregated[family_name] = {
                'cf_samples': cf,
                'weights': w,
                'pdf': pdf_result['pdf'],
                'cdf': pdf_result['cdf'],
                'bin_centers': pdf_result['bin_centers'],
                'bin_edges': pdf_result['bin_edges'],
                'percentiles': percentile_result,
                'statistics': stats,
                'components': [{
                    'pack_code': comp.get('pack_code', 'unknown'),
                    'polymer': comp.get('polymer', 'unknown'),
                    'family': family_name,
                    'n_samples': len(cf)
                }],
                'n_components': 1,
                'n_combined_samples': len(cf)
            }
        else:
            # Multiple components: tensor product aggregation
            # Add family info to each component
            for comp in component_list:
                comp['family'] = family_name

            aggregated[family_name] = aggregate_packaging_components(
                component_list, resolution
            )

    return aggregated
