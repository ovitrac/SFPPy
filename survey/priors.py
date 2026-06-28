"""
survey/priors.py — Triangular Prior Discretization
===================================================

Deterministic discretization of triangular distributions
for survey-scale estimation (no Monte Carlo).

D0 revision (2026-04-08, O. Vitrac):
    - Triangular distribution generalized to (a, c, b) with arbitrary a >= 0,
      replacing the previous hardcoded Triangular(0, mode, max) assumption.
    - Explicit branches for degenerate one-sided triangulars (a = c < b or
      a < c = b).
    - Optional log-spaced low/high segments for wide-span priors (D6 of the
      singleton cp0 revision plan).
    - Node / weight convention: nodes = bin conditional means,
      weights = exact triangular CDF increments — yields exact first moment
      and exact normalisation.

Backward compatibility:
    Any call with min_val = 0.0 (default) reproduces pre-D0 behaviour to
    within the conditional-mean vs midpoint node convention. The overall
    first moment is exact under D0; previously it was approximate.

@project: SFPPy — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

import numpy as np
from typing import Tuple
from survey.models import PriorSpec


# -----------------------------------------------------------------------------
# Triangular CDF and conditional-mean helpers — generalised (a, c, b)
# -----------------------------------------------------------------------------

def triangular_cdf(x: float, min_val: float, mode: float, max_val: float) -> float:
    """
    CDF of Triangular(a, c, b) distribution with a = min_val, c = mode, b = max_val.

    Formula:
        F(x) = 0                                         if x <= a
        F(x) = (x - a)^2 / ((b - a)(c - a))              if a < x <= c
        F(x) = 1 - (b - x)^2 / ((b - a)(b - c))          if c < x < b
        F(x) = 1                                         if x >= b

    Degenerate cases:
        a = c  (right-triangular): only the right branch applies,
               F(x) = 1 - (b - x)^2 / (b - a)^2 on [a, b].
        c = b  (left-triangular):  only the left branch applies,
               F(x) = (x - a)^2 / (b - a)^2 on [a, b].

    Parameters
    ----------
    x : float
        Point at which to evaluate CDF.
    min_val : float
        Lower bound a.
    mode : float
        Mode c, satisfying a <= c <= b.
    max_val : float
        Upper bound b, strictly greater than a.

    Returns
    -------
    float
        CDF value F(x) in [0, 1].
    """
    if x <= min_val:
        return 0.0
    if x >= max_val:
        return 1.0

    a, c, b = min_val, mode, max_val
    span = b - a

    if c == a:
        # Right-triangular: density decreases from peak at a to 0 at b.
        return 1.0 - ((b - x) ** 2) / (span * span)

    if c == b:
        # Left-triangular: density rises from 0 at a to peak at b.
        return ((x - a) ** 2) / (span * span)

    if x <= c:
        return ((x - a) ** 2) / (span * (c - a))
    return 1.0 - ((b - x) ** 2) / (span * (b - c))


def _triangular_bin_conditional_mean(
    x_lo: float, x_hi: float,
    min_val: float, mode: float, max_val: float,
) -> float:
    """
    Closed-form conditional mean E[C | x_lo <= C <= x_hi] of Triangular(a, c, b).

    Computed as the ratio of the incomplete first moment to the CDF mass
    over [x_lo, x_hi]:

        E[C | bin] = integral_{x_lo}^{x_hi} x f(x) dx  /  (F(x_hi) - F(x_lo))

    The PDF is piecewise linear, so the indefinite integral of x f(x) is
    closed-form on each branch (left ramp, right ramp, or either one of the
    two degenerate one-sided cases).

    Handles degenerate bins with zero mass by returning the bin midpoint —
    this is never reached under the grid construction below, which always
    uses bins with strictly positive probability mass.
    """
    a, c, b = min_val, mode, max_val
    span = b - a

    if span <= 0.0:
        return 0.5 * (x_lo + x_hi)

    def _left_antideriv(x: float) -> float:
        # F_left(x) = (x - a)^2 / (span * (c - a))  on [a, c], c > a
        # so f_left(x) = 2 (x - a) / (span * (c - a))
        # integral x f_left dx = (2 / (span * (c - a))) * (x^3/3 - a x^2/2)
        w = span * (c - a)
        return (2.0 / w) * (x ** 3 / 3.0 - 0.5 * a * x * x)

    def _right_antideriv(x: float) -> float:
        # F_right(x) - F_right(c) = ... with f_right(x) = 2 (b - x) / (span * (b - c))
        # integral x f_right dx = (2 / (span * (b - c))) * (b x^2 / 2 - x^3 / 3)
        w = span * (b - c)
        return (2.0 / w) * (0.5 * b * x * x - x ** 3 / 3.0)

    def _right_tri_antideriv(x: float) -> float:
        # a = c: f(x) = 2 (b - x) / span^2 on [a, b]
        return (2.0 / (span * span)) * (0.5 * b * x * x - x ** 3 / 3.0)

    def _left_tri_antideriv(x: float) -> float:
        # c = b: f(x) = 2 (x - a) / span^2 on [a, b]
        return (2.0 / (span * span)) * (x ** 3 / 3.0 - 0.5 * a * x * x)

    # Degenerate: right-triangular (a = c)
    if c == a:
        num = _right_tri_antideriv(x_hi) - _right_tri_antideriv(x_lo)
        denom = triangular_cdf(x_hi, a, c, b) - triangular_cdf(x_lo, a, c, b)
        if denom <= 0.0:
            return 0.5 * (x_lo + x_hi)
        return num / denom

    # Degenerate: left-triangular (c = b)
    if c == b:
        num = _left_tri_antideriv(x_hi) - _left_tri_antideriv(x_lo)
        denom = triangular_cdf(x_hi, a, c, b) - triangular_cdf(x_lo, a, c, b)
        if denom <= 0.0:
            return 0.5 * (x_lo + x_hi)
        return num / denom

    # Ordinary triangular: bin may straddle the mode c
    if x_hi <= c:
        num = _left_antideriv(x_hi) - _left_antideriv(x_lo)
    elif x_lo >= c:
        num = _right_antideriv(x_hi) - _right_antideriv(x_lo)
    else:
        # Straddle: split at c
        num = (_left_antideriv(c) - _left_antideriv(x_lo)) \
              + (_right_antideriv(x_hi) - _right_antideriv(c))

    denom = triangular_cdf(x_hi, a, c, b) - triangular_cdf(x_lo, a, c, b)
    if denom <= 0.0:
        return 0.5 * (x_lo + x_hi)
    return num / denom


# -----------------------------------------------------------------------------
# Segment builders (linear / log spacing)
# -----------------------------------------------------------------------------

def _segment_edges(lo: float, hi: float, n: int, spacing: str) -> np.ndarray:
    """
    Return n+1 edges spanning [lo, hi] with linear or log spacing.

    For log spacing, lo must be > 0.
    """
    if n < 1:
        raise ValueError(f"Segment must have at least 1 bin; got n={n}")
    if hi <= lo:
        raise ValueError(f"Segment hi must exceed lo; got lo={lo}, hi={hi}")
    if spacing == 'linear':
        return np.linspace(lo, hi, n + 1)
    if spacing == 'log':
        if lo <= 0.0:
            raise ValueError(
                f"Log-spaced segment requires lo > 0; got lo={lo}. "
                f"Scenarios with min=0 cannot use log spacing on the low segment."
            )
        return np.exp(np.linspace(np.log(lo), np.log(hi), n + 1))
    raise ValueError(f"Unknown spacing '{spacing}'; must be 'linear' or 'log'")


def _build_segment(
    lo: float, hi: float, n: int, spacing: str,
    min_val: float, mode: float, max_val: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build one segment (low or high) of the grid.

    Returns (nodes, weights) where:
        - edges = n+1 points spanning [lo, hi] with requested spacing
        - nodes = bin conditional means E[C | edge_i <= C <= edge_{i+1}]
        - weights = exact CDF increments F(edge_{i+1}) - F(edge_i)
    """
    edges = _segment_edges(lo, hi, n, spacing)
    nodes = np.empty(n, dtype=float)
    weights = np.empty(n, dtype=float)
    for i in range(n):
        w = (triangular_cdf(edges[i + 1], min_val, mode, max_val)
             - triangular_cdf(edges[i], min_val, mode, max_val))
        weights[i] = w
        if w > 0.0:
            nodes[i] = _triangular_bin_conditional_mean(
                edges[i], edges[i + 1], min_val, mode, max_val
            )
        else:
            # Empty bin — shouldn't occur inside the support, but guard against it.
            nodes[i] = 0.5 * (edges[i] + edges[i + 1])
    return nodes, weights


# -----------------------------------------------------------------------------
# Public grid constructor (D0 + D6)
# -----------------------------------------------------------------------------

def triangular_grid(
    min_val: float,
    mode: float,
    max_val: float,
    n_low: int,
    n_high: int,
    spacing_low: str = 'linear',
    spacing_high: str = 'linear',
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Discretize Triangular(min_val, mode, max_val) into nodes + cell-integrated weights.

    Three branches:
        - min_val == mode  (right-triangular): only the high segment is built,
          n_high nodes covering [min_val, max_val].
        - mode == max_val  (left-triangular):  only the low segment is built,
          n_low nodes covering [min_val, max_val].
        - otherwise: two segments joined at the mode, n_low + n_high nodes.

    Nodes are bin conditional means; weights are exact CDF increments. This
    yields an exact first moment (sum w_i node_i == (a + c + b) / 3) and
    exact normalisation (sum w_i == 1).

    Parameters
    ----------
    min_val : float
        Lower bound a >= 0.
    mode : float
        Mode c, a <= c <= b, c < b unless right-triangular.
    max_val : float
        Upper bound b > a.
    n_low : int
        Number of bins in the low segment [a, c] when a < c.
        Ignored when a == c (right-triangular).
    n_high : int
        Number of bins in the high segment [c, b] when c < b.
        Ignored when c == b (left-triangular).
    spacing_low : str
        'linear' (default) or 'log'. Log requires min_val > 0.
    spacing_high : str
        'linear' (default) or 'log'. Log requires mode > 0.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (nodes, weights), both 1D. Weights sum to 1.0 exactly.
        Length is n_low + n_high for ordinary triangulars; n_high for
        right-triangulars; n_low for left-triangulars.

    Raises
    ------
    ValueError
        If min_val, mode, max_val violate 0 <= a <= c <= b with a < b,
        or if grid parameters are invalid, or if log spacing is requested
        on a segment whose lower bound is zero.
    """
    if not (0.0 <= min_val <= mode <= max_val):
        raise ValueError(
            f"Must have 0 <= min_val <= mode <= max_val. "
            f"Got min_val={min_val}, mode={mode}, max_val={max_val}"
        )
    if min_val == max_val:
        raise ValueError(
            f"Degenerate triangular: min_val == max_val == {min_val}. "
            f"The distribution collapses to a point mass and cannot be discretized."
        )

    # Right-triangular: no low segment
    if min_val == mode:
        if n_high < 1:
            raise ValueError(f"Right-triangular requires n_high >= 1; got n_high={n_high}")
        nodes, weights = _build_segment(
            min_val, max_val, n_high, spacing_high,
            min_val, mode, max_val,
        )
        # Normalise (should be exactly 1 to machine precision, but enforce).
        total = weights.sum()
        if total <= 0.0:
            raise ValueError("Degenerate weights: total mass <= 0.")
        weights /= total
        return nodes, weights

    # Left-triangular: no high segment
    if mode == max_val:
        if n_low < 1:
            raise ValueError(f"Left-triangular requires n_low >= 1; got n_low={n_low}")
        nodes, weights = _build_segment(
            min_val, max_val, n_low, spacing_low,
            min_val, mode, max_val,
        )
        total = weights.sum()
        if total <= 0.0:
            raise ValueError("Degenerate weights: total mass <= 0.")
        weights /= total
        return nodes, weights

    # Ordinary triangular: two segments joined at mode
    if n_low < 1 or n_high < 1:
        raise ValueError(
            f"Ordinary triangular requires n_low >= 1 and n_high >= 1; "
            f"got n_low={n_low}, n_high={n_high}"
        )
    low_nodes, low_w = _build_segment(
        min_val, mode, n_low, spacing_low,
        min_val, mode, max_val,
    )
    high_nodes, high_w = _build_segment(
        mode, max_val, n_high, spacing_high,
        min_val, mode, max_val,
    )
    nodes = np.concatenate([low_nodes, high_nodes])
    weights = np.concatenate([low_w, high_w])
    total = weights.sum()
    if total <= 0.0:
        raise ValueError("Degenerate weights: total mass <= 0.")
    weights /= total
    return nodes, weights


# -----------------------------------------------------------------------------
# Public API bound to PriorSpec
# -----------------------------------------------------------------------------

def discretize_prior(prior: PriorSpec) -> Tuple[np.ndarray, np.ndarray]:
    """
    Discretize a PriorSpec into grid values and weights.

    Honours prior.min_val (D0), prior.spacing_low and prior.spacing_high (D6).

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
        min_val=prior.min_val,
        mode=prior.mode,
        max_val=prior.max_val,
        n_low=prior.n_low,
        n_high=prior.n_high,
        spacing_low=prior.spacing_low,
        spacing_high=prior.spacing_high,
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
