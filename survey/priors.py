"""
survey/priors.py — Triangular Prior Discretization
===================================================

Deterministic discretization of triangular distributions
for survey-scale estimation (no Monte Carlo).

@project: SFPPy/INSERM — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

import numpy as np
from typing import Tuple
from survey.models import PriorSpec


def triangular_cdf(x: float, mode: float, max_val: float) -> float:
    """
    CDF of Triangular(0, mode, max_val) distribution.

    Parameters
    ----------
    x : float
        Point at which to evaluate CDF.
    mode : float
        Mode of triangular distribution.
    max_val : float
        Maximum value.

    Returns
    -------
    float
        CDF value F(x).
    """
    if x <= 0.0:
        return 0.0
    if x >= max_val:
        return 1.0
    if x <= mode:
        return (x * x) / (max_val * mode)
    return 1.0 - ((max_val - x) ** 2) / ((max_val - mode) * max_val)


def triangular_grid(
    mode: float,
    max_val: float,
    n_low: int,
    n_high: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Discretize Triangular(0, mode, max_val) into nodes + cell-integrated weights.

    Creates a grid with n_low points in [0, mode] and n_high points in [mode, max_val].
    Weights are computed as the probability mass in each cell (finite-difference quadrature).

    Parameters
    ----------
    mode : float
        Mode of the triangular distribution.
    max_val : float
        Maximum value of the distribution.
    n_low : int
        Number of grid points in [0, mode].
    n_high : int
        Number of grid points in [mode, max_val].

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (values, weights) where values are grid points and weights sum to 1.

    Raises
    ------
    ValueError
        If mode not in (0, max_val] or grid parameters invalid.
    """
    if not (0.0 < mode <= max_val):
        raise ValueError(f"mode must be in (0, max_val]. Got mode={mode}, max_val={max_val}")
    if n_low < 1 or n_high < 1:
        raise ValueError(f"n_low and n_high must be >= 1. Got n_low={n_low}, n_high={n_high}")

    # Build grid: n_low points from 0 to mode, n_high from mode to max
    low = np.linspace(0.0, mode, num=n_low, endpoint=True)
    high = np.linspace(mode, max_val, num=n_high, endpoint=True)
    values = np.unique(np.concatenate([low, high]))
    values.sort()

    # Build cell edges (midpoints between consecutive values)
    edges = np.empty(values.size + 1)
    edges[0] = 0.0
    edges[-1] = max_val
    for i in range(1, values.size):
        edges[i] = 0.5 * (values[i - 1] + values[i])

    # Compute weights as probability mass in each cell
    weights = np.array([
        triangular_cdf(edges[i + 1], mode, max_val) - triangular_cdf(edges[i], mode, max_val)
        for i in range(values.size)
    ], dtype=float)

    # Normalize (should already be ~1, but ensure numerical stability)
    total = weights.sum()
    if total <= 0:
        raise ValueError("Degenerate triangular weights (sum <= 0)")
    weights /= total

    return values, weights


def discretize_prior(prior: PriorSpec) -> Tuple[np.ndarray, np.ndarray]:
    """
    Discretize a PriorSpec into grid values and weights.

    Parameters
    ----------
    prior : PriorSpec
        Prior specification.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (values, weights).
    """
    return triangular_grid(
        mode=prior.mode,
        max_val=prior.max_val,
        n_low=prior.n_low,
        n_high=prior.n_high,
    )


def compute_weight_matrix(
    time_weights: np.ndarray,
    conc_weights: np.ndarray
) -> np.ndarray:
    """
    Compute joint weight matrix from marginal weights.

    Assumes independence: P(t, c) = P(t) * P(c).

    Parameters
    ----------
    time_weights : np.ndarray
        Weights for time prior (shape: N_T).
    conc_weights : np.ndarray
        Weights for concentration prior (shape: N_C).

    Returns
    -------
    np.ndarray
        Joint weight matrix (shape: N_T x N_C).
    """
    return np.outer(time_weights, conc_weights)
