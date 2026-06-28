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

@project: SFPPy — Survey-scale exposure estimation
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


# =========================================================================
# R5b primitives — shared contact-time aggregation (Phase 2.1)
# =========================================================================
#
# Context: the Phase-1 NPZ contract persists each source's CF as a
# structured tensor aligned with its own axis grids (time, time2, conc).
# For a food × family × simulant group, all sources share the same t1
# (and the same t2 if every source is step-2-active) because the food
# record has one contact event. R5b mandates that the cross-source sum
# must happen ON the shared t axis BEFORE integrating time out — not
# as a tensor product over time.
#
# Primitives below implement:
#   - sources_share_t_axis()          — invariant check
#   - sum_sources_at_fixed_t()        — cross-source sum at a fixed t
#                                       cell (Cp0 tensor product over
#                                       independent per-source priors)
#   - combine_sources_shared_t()      — full pipeline producing the
#                                       weighted empirical CF for a
#                                       food × family × simulant group
#
# Memory note: for K sources each with n_cp0 Cp0 nodes, the per-cell
# Cp0 tensor product has n_cp0^K weighted samples. For K=2 (typical
# bilayer), n_cp0=31: ≈ 1k/cell → ≈ 30k for a 1-step group, ≈ 1M for
# a 2-step group. Tractable. For K=4 (hypothetical B33-like with all
# four layer-sources aggregated directly), n_cp0^4 ≈ 1M/cell ≈ 30M
# samples — tight but workable. Deterministic compression hooks
# (workplan § 6.5 G6) are not implemented here; they are deferred
# until a real dataset triggers a memory breach.


def sources_share_t_axis(sources) -> bool:
    """
    Return True if every source in `sources` has identical `time_vals`,
    and — when any source is step2-active — every step2-active source
    has identical `time2_vals`. Raise ValueError otherwise.

    Each source is a dict with keys:
      time_vals     : np.ndarray
      time_weights  : np.ndarray
      (2-step only) time2_vals, time2_weights : np.ndarray
      step2_active  : Optional[bool]     # True for 2-step sources,
                                          # False / absent otherwise
    """
    if not sources:
        raise ValueError("sources must be non-empty")

    ref_t1 = np.asarray(sources[0]["time_vals"], dtype=float)
    ref_t1_w = np.asarray(sources[0]["time_weights"], dtype=float)
    for i, s in enumerate(sources[1:], start=1):
        t1 = np.asarray(s["time_vals"], dtype=float)
        t1_w = np.asarray(s["time_weights"], dtype=float)
        if t1.shape != ref_t1.shape or not np.allclose(t1, ref_t1):
            raise ValueError(
                f"time_vals mismatch in source {i} "
                f"(shape {t1.shape} vs {ref_t1.shape} or values differ)")
        if t1_w.shape != ref_t1_w.shape or not np.allclose(t1_w, ref_t1_w):
            raise ValueError(f"time_weights mismatch in source {i}")

    active = [s for s in sources if s.get("step2_active")]
    if active:
        ref_t2 = np.asarray(active[0]["time2_vals"], dtype=float)
        ref_t2_w = np.asarray(active[0]["time2_weights"], dtype=float)
        for i, s in enumerate(active[1:], start=1):
            t2 = np.asarray(s["time2_vals"], dtype=float)
            t2_w = np.asarray(s["time2_weights"], dtype=float)
            if t2.shape != ref_t2.shape or not np.allclose(t2, ref_t2):
                raise ValueError(
                    f"time2_vals mismatch in step2-active source {i}")
            if t2_w.shape != ref_t2_w.shape or not np.allclose(t2_w, ref_t2_w):
                raise ValueError(
                    f"time2_weights mismatch in step2-active source {i}")
    return True


def _cf_at_fixed_t(source, t1_idx: int, t2_idx=None) -> np.ndarray:
    """
    Slice the source's CF_tensor at a fixed (t1_idx, t2_idx).

    Returns a 1-D array of CF values indexed by the source's Cp0 grid.

    For a step-2-active source with 3-D CF_tensor (n_t1, n_t2, n_cp0),
    this returns CF_tensor[t1_idx, t2_idx, :].

    For a step-2-inactive source (scenario has no step 2 OR this source
    is a t2=0 "removed" component), CF_tensor is 2-D (n_t1, n_cp0) and
    t2_idx is ignored — the source contributes the same CF regardless
    of t2, which is the correct physics (no step 2 for this source).
    """
    cf = np.asarray(source["CF_tensor"])
    if cf.ndim == 3:
        if t2_idx is None:
            raise ValueError(
                "t2_idx required for step-2-active source with 3-D CF_tensor")
        return cf[t1_idx, t2_idx, :]
    elif cf.ndim == 2:
        return cf[t1_idx, :]
    raise ValueError(f"unsupported CF_tensor.ndim = {cf.ndim}")


def sum_sources_at_fixed_t(sources, t1_idx: int, t2_idx=None):
    """
    Cross-source CF sum at a fixed (t1_idx, t2_idx) cell.

    Each source contributes a 1-D CF(Cp0) vector (sliced from its
    CF_tensor at the given t indices) and an independent Cp0 prior.
    The sum-of-independent-random-variables is the tensor product of
    (CF_s, w_cp0_s) pairs — `combine_multiple_tensors` is exactly this.

    Returns
    -------
    (cf_samples, weights) : Tuple[np.ndarray, np.ndarray]
        Flattened weighted empirical distribution of the cross-source
        CF sum at the fixed t cell, shape (Π n_cp0_s,).

    Notes
    -----
    Step-2-active and step-2-inactive sources can be mixed freely in
    `sources`. Inactive sources' 2-D CF_tensor is sliced with t1_idx
    only, so their contribution is constant across t2 — which is the
    correct physical semantics of a removed-at-step-2 component.
    """
    tensors = []
    for s in sources:
        cf_slice = _cf_at_fixed_t(s, t1_idx, t2_idx)
        w_cp0 = np.asarray(s["conc_weights"], dtype=float)
        tensors.append((cf_slice, w_cp0))
    return combine_multiple_tensors(tensors)


# Above this exact flat-sample count the group is aggregated by Monte-Carlo
# instead of the exact tensor product. The exact path enumerates
#   n_t1 * (n_t2 or 1) * prod_s(n_cp0_s)
# samples — prod_s(n_cp0_s) = 30^K for K sources, so it explodes for K>=3
# (generic family_ids like 'monomer' span many pack components → large K).
MAX_EXACT_FLAT_SAMPLES = 2_000_000
MC_SAMPLES = 200_000
MC_SEED = 20260523  # deterministic → reproducible aggregation


def _norm_p(w: np.ndarray) -> np.ndarray:
    w = np.asarray(w, dtype=float).ravel()
    w = np.clip(w, 0.0, None)
    s = w.sum()
    return (w / s) if s > 0 else np.full(len(w), 1.0 / max(len(w), 1))


def _exact_flat_size(sources) -> float:
    """Flat sample count the exact tensor-product path would materialise."""
    ref = sources[0]
    n_t1 = len(np.asarray(ref["time_weights"]).ravel())
    active = [s for s in sources if s.get("step2_active")]
    n_t2 = len(np.asarray(active[0]["time2_weights"]).ravel()) if active else 1
    prod_cp0 = 1.0
    for s in sources:
        prod_cp0 *= len(np.asarray(s["conc_weights"]).ravel())
    return n_t1 * n_t2 * prod_cp0


def combine_sources_shared_t_mc(sources, n_mc: int = MC_SAMPLES,
                                seed: int = MC_SEED):
    """Monte-Carlo equivalent of `combine_sources_shared_t` for large groups.

    Draws `n_mc` joint samples honoring R5b shared contact time: each sample
    draws ONE (t1, t2) index from the shared time weights (applied to every
    source), and an INDEPENDENT Cp0 index per source, then sums the per-source
    CF slices. O(n_mc · K) and O(n_mc) memory — scales to any group size.

    Deterministic (fixed seed) so the aggregation is reproducible. Returns
    (cf_samples, weights) with uniform weights, ready for the percentile/
    statistics helpers exactly like the exact path.
    """
    if not sources:
        raise ValueError("sources must be non-empty")
    sources_share_t_axis(sources)
    rng = np.random.default_rng(seed)

    ref = sources[0]
    t1_p = _norm_p(ref["time_weights"])
    t1_idx = rng.choice(len(t1_p), size=n_mc, p=t1_p)

    active = [s for s in sources if s.get("step2_active")]
    if active:
        t2_p = _norm_p(active[0]["time2_weights"])
        t2_idx = rng.choice(len(t2_p), size=n_mc, p=t2_p)
    else:
        t2_idx = None

    cf_total = np.zeros(n_mc, dtype=float)
    for s in sources:
        cf = np.asarray(s["CF_tensor"], dtype=float)
        cp0_p = _norm_p(s["conc_weights"])
        cp0_idx = rng.choice(len(cp0_p), size=n_mc, p=cp0_p)
        if cf.ndim == 3:
            if t2_idx is None:
                raise ValueError("3-D source in a group with no step-2 time axis")
            cf_total += cf[t1_idx, t2_idx, cp0_idx]
        elif cf.ndim == 2:
            cf_total += cf[t1_idx, cp0_idx]   # step-2-inactive: constant in t2
        else:
            raise ValueError(f"unsupported CF_tensor.ndim = {cf.ndim}")

    return cf_total, np.full(n_mc, 1.0 / n_mc, dtype=float)


def combine_sources_shared_t(sources):
    """
    Aggregate all sources of a food × family × simulant group into a
    single weighted empirical CF distribution, honoring R5b shared
    contact time.

    Dispatch: groups whose exact tensor-product would stay within
    ``MAX_EXACT_FLAT_SAMPLES`` use the exact enumeration below (unchanged,
    preserving prior results); larger groups fall back to the deterministic
    Monte-Carlo path ``combine_sources_shared_t_mc`` so aggregation scales to
    142k+ jobs without the 30^K blow-up.

    Pipeline:
      1. verify shared t axes (raise if violated)
      2. for each (t1_i, t2_j) cell on the shared grid:
           cf_cell, w_cell = sum_sources_at_fixed_t(...)
      3. accumulate with shared time weights:
           weighted samples = cf_cell, w_cell * w_t1[i] * w_t2[j]
      4. flatten and return as (cf_samples, weights)

    The returned distribution is normalised (sum of weights = 1) so
    it can be fed directly to compute_percentiles / compute_risk_
    percentiles / compute_statistics.

    Parameters
    ----------
    sources : List[dict]
        Each dict has keys: CF_tensor, time_vals, time_weights,
        conc_vals, conc_weights, and optionally time2_vals,
        time2_weights, step2_active.

    Returns
    -------
    (cf_samples, weights) : Tuple[np.ndarray, np.ndarray]
        Flattened weighted empirical CF distribution for the group.
    """
    if not sources:
        raise ValueError("sources must be non-empty")

    sources_share_t_axis(sources)

    # Scale guard: fall back to Monte-Carlo for groups that would explode the
    # exact tensor product (K >= 3 sources). Small groups stay exact.
    if _exact_flat_size(sources) > MAX_EXACT_FLAT_SAMPLES:
        return combine_sources_shared_t_mc(sources)

    ref = sources[0]
    t1_w = np.asarray(ref["time_weights"], dtype=float)
    n_t1 = len(t1_w)

    active = [s for s in sources if s.get("step2_active")]
    has_step2 = bool(active)
    if has_step2:
        t2_w = np.asarray(active[0]["time2_weights"], dtype=float)
        n_t2 = len(t2_w)
    else:
        t2_w = None
        n_t2 = None

    cf_chunks: List[np.ndarray] = []
    w_chunks: List[np.ndarray] = []

    if has_step2:
        for i in range(n_t1):
            for j in range(n_t2):
                cf_ij, w_ij = sum_sources_at_fixed_t(sources, i, j)
                cf_chunks.append(np.asarray(cf_ij, dtype=float))
                w_chunks.append(np.asarray(w_ij, dtype=float)
                                * t1_w[i] * t2_w[j])
    else:
        for i in range(n_t1):
            cf_i, w_i = sum_sources_at_fixed_t(sources, i)
            cf_chunks.append(np.asarray(cf_i, dtype=float))
            w_chunks.append(np.asarray(w_i, dtype=float) * t1_w[i])

    cf_flat = np.concatenate(cf_chunks)
    w_flat = np.concatenate(w_chunks)

    total = w_flat.sum()
    if total > 0:
        w_flat = w_flat / total
    return cf_flat, w_flat


def jaccard_similarity(cas_sets) -> float:
    """
    Minimum pairwise Jaccard similarity over a collection of CAS sets.

    Returns 1.0 for singletons or if all sets are identical;
    returns a value in [0, 1] reflecting the worst overlap across
    all pairs of sources. Used as a diagnostic flag per R5b (§ 6.2 of
    multicomponent_bilayer_aggregation_20260415.md) — it does NOT
    gate summation.
    """
    cas_sets = [set(cs) for cs in cas_sets]
    if len(cas_sets) <= 1:
        return 1.0
    j_min = 1.0
    for i in range(len(cas_sets)):
        for k in range(i + 1, len(cas_sets)):
            a, b = cas_sets[i], cas_sets[k]
            union = a | b
            inter = a & b
            j = (len(inter) / len(union)) if union else 1.0
            if j < j_min:
                j_min = j
    return j_min
