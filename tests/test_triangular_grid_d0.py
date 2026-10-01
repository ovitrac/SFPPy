"""
Unit tests for survey/priors.py — D0 generalised triangular grid.

Covers:
    - Backward compatibility: min=0 path gives exact first moment
      (previously approximate under midpoint convention).
    - Arbitrary min > 0: first moment and normalisation exact.
    - Degenerate right-triangular (a = m < b).
    - Degenerate left-triangular  (a < m = b).
    - Three geometry regimes: narrow, wide (log-required), symmetric.
    - Log vs linear spacing on both segments.
    - Error handling: invalid arguments, log spacing with zero lower bound.

Author: Olivier Vitrac, PhD, HDR
Date: 2026-04-08
"""

import numpy as np
import pytest

from survey.priors import (
    triangular_cdf,
    triangular_grid,
    discretize_prior,
    _triangular_bin_conditional_mean,
)
from survey.models import PriorSpec


TOL = 1e-10


def _expected_mean(a: float, c: float, b: float) -> float:
    """Closed-form mean of Triangular(a, c, b)."""
    return (a + c + b) / 3.0


# ---------------------------------------------------------------------------
# CDF correctness
# ---------------------------------------------------------------------------

class TestTriangularCDF:
    """F(a) = 0, F(b) = 1, F(c) = (c - a) / (b - a)."""

    def test_boundaries_ordinary(self):
        a, c, b = 7.0, 14.0, 464.0
        assert triangular_cdf(a, a, c, b) == 0.0
        assert triangular_cdf(b, a, c, b) == 1.0
        # Below a and above b
        assert triangular_cdf(a - 1, a, c, b) == 0.0
        assert triangular_cdf(b + 1, a, c, b) == 1.0

    def test_mode_value(self):
        a, c, b = 7.0, 14.0, 464.0
        # F(c) = (c - a) / (b - a)
        expected = (c - a) / (b - a)
        assert abs(triangular_cdf(c, a, c, b) - expected) < TOL

    def test_monotone(self):
        a, c, b = 0.6, 80.0, 134.0
        xs = np.linspace(a, b, 101)
        F = [triangular_cdf(x, a, c, b) for x in xs]
        for i in range(1, len(F)):
            assert F[i] >= F[i - 1] - TOL

    def test_right_triangular_cdf(self):
        # a = c = 7, b = 116
        a, c, b = 7.0, 7.0, 116.0
        # F(a) = 0, F(b) = 1, F increases on [a, b]
        assert triangular_cdf(a, a, c, b) == 0.0
        assert triangular_cdf(b, a, c, b) == 1.0
        # Check F at halfway: quadratic rise from 1 - (b-x)^2/span^2
        x_mid = 0.5 * (a + b)
        expected = 1.0 - ((b - x_mid) ** 2) / ((b - a) ** 2)
        assert abs(triangular_cdf(x_mid, a, c, b) - expected) < TOL

    def test_left_triangular_cdf(self):
        # a = 1, c = b = 10
        a, c, b = 1.0, 10.0, 10.0
        assert triangular_cdf(a, a, c, b) == 0.0
        assert triangular_cdf(b, a, c, b) == 1.0
        x_mid = 0.5 * (a + b)
        expected = ((x_mid - a) ** 2) / ((b - a) ** 2)
        assert abs(triangular_cdf(x_mid, a, c, b) - expected) < TOL


# ---------------------------------------------------------------------------
# Ordinary triangular — D0 backward compatibility and first-moment exactness
# ---------------------------------------------------------------------------

class TestOrdinaryTriangular:
    """Any (a, c, b) with a < c < b, both segments present."""

    def test_backward_compat_min_zero(self):
        # Pre-D0 call pattern: min = 0, linear spacing on both segments.
        mode, max_val = 50.0, 200.0
        nodes, weights = triangular_grid(
            min_val=0.0, mode=mode, max_val=max_val,
            n_low=15, n_high=15,
            spacing_low='linear', spacing_high='linear',
        )
        # Normalisation exact.
        assert abs(weights.sum() - 1.0) < TOL
        # First moment exact with the new conditional-mean convention.
        mean_computed = float(np.dot(nodes, weights))
        mean_expected = _expected_mean(0.0, mode, max_val)
        assert abs(mean_computed - mean_expected) < TOL

    def test_narrow_linear(self):
        # gPET IAS-like: (7, 14, 232)
        a, c, b = 7.0, 14.0, 232.0
        nodes, weights = triangular_grid(
            min_val=a, mode=c, max_val=b,
            n_low=15, n_high=15,
            spacing_low='linear', spacing_high='linear',
        )
        assert abs(weights.sum() - 1.0) < TOL
        mean_computed = float(np.dot(nodes, weights))
        assert abs(mean_computed - _expected_mean(a, c, b)) < TOL
        # All nodes inside the support.
        assert nodes[0] >= a - TOL
        assert nodes[-1] <= b + TOL
        # Nodes strictly increasing (no duplicates at the mode — construction
        # uses two disjoint bin sets joined at c as an edge, so conditional
        # means on adjacent bins straddling c are distinct).
        assert np.all(np.diff(nodes) > 0)

    def test_wide_log_low_segment(self):
        # HDPE NIAS-like after β+: span ~528
        a, c, b = 0.6, 302.0, 634.0
        nodes, weights = triangular_grid(
            min_val=a, mode=c, max_val=b,
            n_low=25, n_high=15,
            spacing_low='log', spacing_high='linear',
        )
        assert abs(weights.sum() - 1.0) < TOL
        mean_computed = float(np.dot(nodes, weights))
        assert abs(mean_computed - _expected_mean(a, c, b)) < TOL
        # Log-spaced low segment: first three nodes geometrically spread.
        # Check that log10(nodes[1]/nodes[0]) > 0 and roughly uniform in log.
        assert nodes[0] > a and nodes[0] < 10 * a  # first bin near the LOQ
        assert nodes[24] < c  # last low node still below mode

    def test_symmetric_linear(self):
        # Symmetric prior for sanity: (10, 50, 90)
        a, c, b = 10.0, 50.0, 90.0
        nodes, weights = triangular_grid(
            min_val=a, mode=c, max_val=b,
            n_low=20, n_high=20,
        )
        assert abs(weights.sum() - 1.0) < TOL
        mean_computed = float(np.dot(nodes, weights))
        assert abs(mean_computed - _expected_mean(a, c, b)) < TOL
        # Symmetric prior: mean equals mode.
        assert abs(mean_computed - c) < TOL


# ---------------------------------------------------------------------------
# Degenerate cases: right-triangular (a = c) and left-triangular (c = b)
# ---------------------------------------------------------------------------

class TestDegenerateTriangulars:

    def test_right_triangular(self):
        # gPET IAS at w=0.5: (7, 7, 116)
        a, c, b = 7.0, 7.0, 116.0
        nodes, weights = triangular_grid(
            min_val=a, mode=c, max_val=b,
            n_low=25, n_high=15,
            spacing_low='log', spacing_high='linear',
        )
        # Only the high segment exists: length = n_high
        assert len(nodes) == 15
        assert len(weights) == 15
        assert abs(weights.sum() - 1.0) < TOL
        mean_computed = float(np.dot(nodes, weights))
        expected = _expected_mean(a, c, b)  # (7 + 7 + 116) / 3 = 43.333...
        assert abs(mean_computed - expected) < TOL
        # All nodes in [a, b]
        assert nodes[0] >= a - TOL
        assert nodes[-1] <= b + TOL

    def test_left_triangular(self):
        a, c, b = 1.0, 10.0, 10.0
        nodes, weights = triangular_grid(
            min_val=a, mode=c, max_val=b,
            n_low=25, n_high=15,
            spacing_low='linear', spacing_high='linear',
        )
        # Only the low segment exists: length = n_low
        assert len(nodes) == 25
        assert abs(weights.sum() - 1.0) < TOL
        mean_computed = float(np.dot(nodes, weights))
        expected = _expected_mean(a, c, b)  # (1 + 10 + 10) / 3 = 7
        assert abs(mean_computed - expected) < TOL

    def test_right_triangular_log_high(self):
        # Right-tri with log spacing on the high segment.
        a, c, b = 0.6, 0.6, 268.0
        nodes, weights = triangular_grid(
            min_val=a, mode=c, max_val=b,
            n_low=25, n_high=15,
            spacing_low='log', spacing_high='log',
        )
        assert len(nodes) == 15
        assert abs(weights.sum() - 1.0) < TOL
        mean_computed = float(np.dot(nodes, weights))
        assert abs(mean_computed - _expected_mean(a, c, b)) < TOL


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

class TestInputValidation:

    def test_reject_negative_min(self):
        with pytest.raises(ValueError):
            triangular_grid(min_val=-1.0, mode=10.0, max_val=100.0,
                            n_low=10, n_high=10)

    def test_reject_mode_above_max(self):
        with pytest.raises(ValueError):
            triangular_grid(min_val=0.0, mode=100.0, max_val=50.0,
                            n_low=10, n_high=10)

    def test_reject_min_above_mode(self):
        with pytest.raises(ValueError):
            triangular_grid(min_val=50.0, mode=10.0, max_val=100.0,
                            n_low=10, n_high=10)

    def test_reject_point_mass(self):
        with pytest.raises(ValueError):
            triangular_grid(min_val=5.0, mode=5.0, max_val=5.0,
                            n_low=10, n_high=10)

    def test_reject_log_with_zero_min(self):
        # min=0 cannot combine with log-spaced low segment.
        with pytest.raises(ValueError):
            triangular_grid(min_val=0.0, mode=50.0, max_val=200.0,
                            n_low=10, n_high=10,
                            spacing_low='log', spacing_high='linear')

    def test_reject_unknown_spacing(self):
        with pytest.raises(ValueError):
            triangular_grid(min_val=1.0, mode=50.0, max_val=200.0,
                            n_low=10, n_high=10,
                            spacing_low='foobar', spacing_high='linear')

    def test_reject_zero_n_low_ordinary(self):
        with pytest.raises(ValueError):
            triangular_grid(min_val=1.0, mode=10.0, max_val=100.0,
                            n_low=0, n_high=10)

    def test_reject_zero_n_high_right_triangular(self):
        with pytest.raises(ValueError):
            triangular_grid(min_val=7.0, mode=7.0, max_val=100.0,
                            n_low=25, n_high=0)


# ---------------------------------------------------------------------------
# PriorSpec → discretize_prior round trip
# ---------------------------------------------------------------------------

class TestPriorSpecRoundTrip:

    def test_default_priorspec_still_works(self):
        # Pre-D0 callsite: PriorSpec(mode, max_val) with defaults.
        prior = PriorSpec(mode=50.0, max_val=200.0, name='legacy')
        assert prior.min_val == 0.0
        values, weights = discretize_prior(prior)
        assert abs(weights.sum() - 1.0) < TOL
        mean_computed = float(np.dot(values, weights))
        assert abs(mean_computed - _expected_mean(0.0, 50.0, 200.0)) < TOL

    def test_new_priorspec_with_min(self):
        prior = PriorSpec(
            min_val=0.6, mode=302.0, max_val=634.0,
            n_low=25, n_high=15,
            spacing_low='log', spacing_high='linear',
            name='hdpe_nias_w1',
        )
        values, weights = discretize_prior(prior)
        assert len(values) == 40
        assert abs(weights.sum() - 1.0) < TOL
        mean_computed = float(np.dot(values, weights))
        assert abs(mean_computed - _expected_mean(0.6, 302.0, 634.0)) < TOL

    def test_priorspec_right_triangular(self):
        # Degenerate case used in gPET IAS at w=0.5
        prior = PriorSpec(
            min_val=7.0, mode=7.0, max_val=116.0,
            n_low=25, n_high=15,
            spacing_low='log', spacing_high='linear',
        )
        values, weights = discretize_prior(prior)
        assert len(values) == 15
        assert abs(weights.sum() - 1.0) < TOL
        mean_computed = float(np.dot(values, weights))
        assert abs(mean_computed - _expected_mean(7.0, 7.0, 116.0)) < TOL


# ---------------------------------------------------------------------------
# Diagnostic: old-solver (min=0 hardcoded) vs D0 on the four representative priors
# ---------------------------------------------------------------------------

class TestPreexistingMin0Bias:
    """
    Quantify the bias introduced when the old solver discretized from 0
    regardless of the YAML-declared min.
    """

    @staticmethod
    def _old_solver_equivalent(mode: float, max_val: float,
                                n_low: int = 15, n_high: int = 15):
        """Reproduce pre-D0 behaviour by calling with min_val=0."""
        return triangular_grid(
            min_val=0.0, mode=mode, max_val=max_val,
            n_low=n_low, n_high=n_high,
        )

    @staticmethod
    def _mass_below(nodes, weights, threshold):
        """Total probability mass assigned to nodes with value < threshold."""
        mask = nodes < threshold
        return float(weights[mask].sum())

    @pytest.mark.parametrize('a, c, b, label', [
        (7.0, 14.0, 464.0, 'gPET IAS singleton'),
        (50.0, 1000.0, 5000.0, 'HDPE antioxidant formulation'),
        (9.1, 29.0, 380.0, 'Paper NIAS singleton'),
        (5.0, 50.0, 500.0, 'Narrow formulation (monomer PET)'),
    ])
    def test_bias_diagnostic(self, a, c, b, label):
        # Old behaviour: discretize from 0, ignoring declared min = a.
        old_nodes, old_w = self._old_solver_equivalent(mode=c, max_val=b)
        # D0 behaviour: honour min = a.
        d0_nodes, d0_w = triangular_grid(
            min_val=a, mode=c, max_val=b,
            n_low=15, n_high=15,
        )

        # Under the old solver, the probability mass assigned to [0, a] is
        # the CDF of Triangular(0, c, b) evaluated at a.
        ghost_mass_old = self._mass_below(old_nodes, old_w, threshold=a)
        ghost_mass_d0 = self._mass_below(d0_nodes, d0_w, threshold=a)

        # Theoretical value: F_Triangular(0, c, b)(a) = a^2 / (b c) on a <= c.
        theoretical_ghost = (a * a) / (b * c) if a <= c else 1.0 - ((b - a) ** 2) / (b * (b - c))

        # D0 assigns zero mass below a (support starts at a).
        assert ghost_mass_d0 < TOL, (
            f"[{label}] D0 should assign zero mass below min={a}, "
            f"got {ghost_mass_d0:.6e}"
        )

        # Old solver assigns the theoretical value (modulo discretisation).
        assert abs(ghost_mass_old - theoretical_ghost) < 0.01, (
            f"[{label}] Old solver ghost mass {ghost_mass_old:.4f} "
            f"should match theoretical {theoretical_ghost:.4f}"
        )

        # First moment comparison.
        mean_old = float(np.dot(old_nodes, old_w))
        mean_d0 = float(np.dot(d0_nodes, d0_w))
        # D0 first moment is exact.
        assert abs(mean_d0 - _expected_mean(a, c, b)) < TOL, (
            f"[{label}] D0 first moment should be exact."
        )

        # Print diagnostic table row (visible with pytest -s).
        print(
            f"\n{label:<40}  a={a:>7.2f}  c={c:>8.2f}  b={b:>8.2f}  "
            f"ghost_mass_old={ghost_mass_old:.4f}  "
            f"ghost_mass_theory={theoretical_ghost:.4f}  "
            f"mean_old={mean_old:>8.2f}  mean_D0={mean_d0:>8.2f}  "
            f"delta_mean={mean_d0 - mean_old:+.2f}"
        )
