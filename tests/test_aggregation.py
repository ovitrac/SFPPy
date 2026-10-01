"""
Unit tests for survey/aggregation.py — Tensor Product Aggregation

Tests the deterministic tensor product approach for multi-component
packaging aggregation.

@project: SFPPy — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@license: MIT
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_array_less


class TestCombineTensors:
    """Tests for combine_tensors() function."""

    def test_basic_combination(self):
        """Two simple tensors should combine correctly."""
        from survey.aggregation import combine_tensors

        cf1 = np.array([1.0, 2.0, 3.0])
        w1 = np.array([0.2, 0.5, 0.3])
        cf2 = np.array([0.5, 1.5])
        w2 = np.array([0.4, 0.6])

        cf_combined, w_combined = combine_tensors(cf1, w1, cf2, w2)

        # Check output shapes
        assert cf_combined.shape == (6,)  # 3 × 2
        assert w_combined.shape == (6,)

        # Check weights sum to 1
        assert_allclose(w_combined.sum(), 1.0, rtol=1e-10)

    def test_cf_values_are_sums(self):
        """Combined CF values should be sums of input values."""
        from survey.aggregation import combine_tensors

        cf1 = np.array([1.0, 2.0])
        w1 = np.array([0.5, 0.5])
        cf2 = np.array([0.1, 0.2])
        w2 = np.array([0.5, 0.5])

        cf_combined, _ = combine_tensors(cf1, w1, cf2, w2)

        # Expected: [1+0.1, 1+0.2, 2+0.1, 2+0.2] = [1.1, 1.2, 2.1, 2.2]
        expected = np.array([1.1, 1.2, 2.1, 2.2])
        assert_allclose(sorted(cf_combined), sorted(expected), rtol=1e-10)

    def test_weights_are_products(self):
        """Combined weights should be products of input weights."""
        from survey.aggregation import combine_tensors

        cf1 = np.array([1.0, 2.0])
        w1 = np.array([0.3, 0.7])
        cf2 = np.array([0.1, 0.2])
        w2 = np.array([0.4, 0.6])

        _, w_combined = combine_tensors(cf1, w1, cf2, w2)

        # Expected weights (before renorm): [0.3×0.4, 0.3×0.6, 0.7×0.4, 0.7×0.6]
        expected_raw = np.array([0.12, 0.18, 0.28, 0.42])
        expected = expected_raw / expected_raw.sum()

        assert_allclose(sorted(w_combined), sorted(expected), rtol=1e-10)

    def test_single_element_tensors(self):
        """Single-element tensors should work correctly."""
        from survey.aggregation import combine_tensors

        cf1 = np.array([5.0])
        w1 = np.array([1.0])
        cf2 = np.array([3.0])
        w2 = np.array([1.0])

        cf_combined, w_combined = combine_tensors(cf1, w1, cf2, w2)

        assert_allclose(cf_combined, [8.0], rtol=1e-10)
        assert_allclose(w_combined, [1.0], rtol=1e-10)

    def test_large_tensors(self):
        """Large tensors (15×15 grid) should handle correctly."""
        from survey.aggregation import combine_tensors

        n = 225  # 15×15 grid
        cf1 = np.linspace(0, 10, n)
        w1 = np.random.rand(n)
        w1 /= w1.sum()

        cf2 = np.linspace(0, 5, n)
        w2 = np.random.rand(n)
        w2 /= w2.sum()

        cf_combined, w_combined = combine_tensors(cf1, w1, cf2, w2)

        assert cf_combined.shape == (n * n,)  # 225 × 225 = 50625
        assert_allclose(w_combined.sum(), 1.0, rtol=1e-10)

        # Max combined should be max(cf1) + max(cf2)
        assert_allclose(cf_combined.max(), cf1.max() + cf2.max(), rtol=1e-10)


class TestCombineMultipleTensors:
    """Tests for combine_multiple_tensors() function."""

    def test_single_tensor(self):
        """Single tensor should be returned unchanged."""
        from survey.aggregation import combine_multiple_tensors

        cf = np.array([1.0, 2.0, 3.0])
        w = np.array([0.2, 0.5, 0.3])

        cf_out, w_out = combine_multiple_tensors([(cf, w)])

        assert_allclose(cf_out, cf, rtol=1e-10)
        assert_allclose(w_out, w, rtol=1e-10)

    def test_three_tensors(self):
        """Three tensors should combine via chained products."""
        from survey.aggregation import combine_multiple_tensors

        cf1 = np.array([1.0, 2.0])
        w1 = np.array([0.5, 0.5])
        cf2 = np.array([0.1])
        w2 = np.array([1.0])
        cf3 = np.array([0.01, 0.02])
        w3 = np.array([0.5, 0.5])

        cf_out, w_out = combine_multiple_tensors([
            (cf1, w1), (cf2, w2), (cf3, w3)
        ])

        # 2 × 1 × 2 = 4 outputs
        assert cf_out.shape == (4,)
        assert_allclose(w_out.sum(), 1.0, rtol=1e-10)

    def test_empty_list_raises(self):
        """Empty list should raise ValueError."""
        from survey.aggregation import combine_multiple_tensors

        with pytest.raises(ValueError):
            combine_multiple_tensors([])


class TestAggregateComponents:
    """Tests for aggregate_components() function."""

    def test_single_component(self):
        """Single component should pass through unchanged."""
        from survey.aggregation import aggregate_components

        component_results = {
            'S1': {
                'CF_samples': np.array([1.0, 2.0, 3.0]),
                'weights': np.array([0.2, 0.5, 0.3]),
                'substance_ids': ['A', 'B'],
            }
        }

        agg = aggregate_components(component_results)

        # Each substance should appear
        assert 'A' in agg
        assert 'B' in agg

        # For single component, values should be unchanged
        assert_allclose(agg['A']['CF_samples'], [1.0, 2.0, 3.0], rtol=1e-10)

    def test_two_components_shared_substance(self):
        """Shared substance should be aggregated via tensor product."""
        from survey.aggregation import aggregate_components

        component_results = {
            'S1': {
                'CF_samples': np.array([1.0, 2.0]),
                'weights': np.array([0.5, 0.5]),
                'substance_ids': ['A'],
            },
            'S2': {
                'CF_samples': np.array([0.1, 0.2]),
                'weights': np.array([0.5, 0.5]),
                'substance_ids': ['A'],
            }
        }

        agg = aggregate_components(component_results)

        # Substance A should have combined results
        assert 'A' in agg
        assert len(agg['A']['source_components']) == 2
        assert 'S1' in agg['A']['source_components']
        assert 'S2' in agg['A']['source_components']

        # Combined should have 2×2 = 4 samples
        assert agg['A']['CF_samples'].shape == (4,)

    def test_unique_substances(self):
        """Substances unique to one component should not be combined."""
        from survey.aggregation import aggregate_components

        component_results = {
            'S1': {
                'CF_samples': np.array([1.0, 2.0]),
                'weights': np.array([0.5, 0.5]),
                'substance_ids': ['A'],
            },
            'S2': {
                'CF_samples': np.array([0.1, 0.2]),
                'weights': np.array([0.5, 0.5]),
                'substance_ids': ['B'],
            }
        }

        agg = aggregate_components(component_results)

        # Both substances should appear
        assert 'A' in agg
        assert 'B' in agg

        # Each should have only 2 samples (no combination)
        assert agg['A']['CF_samples'].shape == (2,)
        assert agg['B']['CF_samples'].shape == (2,)

        # Each should have only one source component
        assert len(agg['A']['source_components']) == 1
        assert len(agg['B']['source_components']) == 1


class TestAggregateFamilyWeighted:
    """Tests for aggregate_family_weighted() function."""

    def test_equal_weights(self):
        """Equal weights should produce uniform mixture."""
        from survey.aggregation import aggregate_family_weighted

        substance_results = {
            'A': {
                'CF_samples': np.array([1.0, 2.0]),
                'weights': np.array([0.5, 0.5]),
            },
            'B': {
                'CF_samples': np.array([3.0, 4.0]),
                'weights': np.array([0.5, 0.5]),
            }
        }
        substance_weights = {'A': 1.0, 'B': 1.0}

        cf, w = aggregate_family_weighted(substance_results, substance_weights)

        # Should have 4 samples total (2 + 2)
        assert cf.shape == (4,)
        assert_allclose(w.sum(), 1.0, rtol=1e-10)

    def test_unequal_weights(self):
        """Unequal occurrence weights should affect relative contributions."""
        from survey.aggregation import aggregate_family_weighted

        substance_results = {
            'A': {
                'CF_samples': np.array([1.0, 2.0]),
                'weights': np.array([0.5, 0.5]),
            },
            'B': {
                'CF_samples': np.array([3.0, 4.0]),
                'weights': np.array([0.5, 0.5]),
            }
        }
        # B has double the occurrence weight
        substance_weights = {'A': 1.0, 'B': 2.0}

        cf, w = aggregate_family_weighted(substance_results, substance_weights)

        # Total weight from B should be ~2/3, from A ~1/3
        b_indices = [i for i, v in enumerate(cf) if v >= 3.0]
        b_weight = sum(w[i] for i in b_indices)

        assert_allclose(b_weight, 2.0 / 3.0, rtol=0.1)

    def test_empty_raises(self):
        """Empty input should raise ValueError."""
        from survey.aggregation import aggregate_family_weighted

        with pytest.raises(ValueError):
            aggregate_family_weighted({}, {})


class TestComputePdfFromSamples:
    """Tests for compute_pdf_from_samples() function."""

    def test_basic_pdf(self):
        """PDF should integrate to approximately 1."""
        from survey.aggregation import compute_pdf_from_samples

        cf = np.linspace(0, 10, 100)
        w = np.ones(100) / 100

        result = compute_pdf_from_samples(cf, w, n_bins=50)

        # PDF should exist and have correct shape
        assert 'pdf' in result
        assert result['pdf'].shape == (50,)

        # CDF should end at ~1
        assert_allclose(result['cdf'][-1], 1.0, atol=0.1)

    def test_all_zeros_edge_case(self):
        """All-zero samples should not crash."""
        from survey.aggregation import compute_pdf_from_samples

        cf = np.zeros(100)
        w = np.ones(100) / 100

        result = compute_pdf_from_samples(cf, w)

        assert 'pdf' in result
        assert 'cdf' in result


class TestQuantileFromSamples:
    """Tests for quantile_from_samples() function."""

    def test_median(self):
        """Median of uniform samples should be middle value."""
        from survey.aggregation import quantile_from_samples

        cf = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        w = np.array([0.2, 0.2, 0.2, 0.2, 0.2])

        q50 = quantile_from_samples(cf, w, 0.5)

        # Should be close to 3.0
        assert_allclose(q50, 3.0, atol=0.5)

    def test_extremes(self):
        """Q0 should be min, Q1 should be max."""
        from survey.aggregation import quantile_from_samples

        cf = np.array([1.0, 5.0, 10.0])
        w = np.array([0.33, 0.34, 0.33])

        q0 = quantile_from_samples(cf, w, 0.0)
        q1 = quantile_from_samples(cf, w, 1.0)

        assert q0 == 1.0
        assert q1 == 10.0


class TestComputeStatistics:
    """Tests for compute_statistics() function."""

    def test_basic_statistics(self):
        """Statistics should be computed correctly."""
        from survey.aggregation import compute_statistics

        cf = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        w = np.array([0.2, 0.2, 0.2, 0.2, 0.2])

        stats = compute_statistics(cf, w)

        assert 'mean' in stats
        assert 'std' in stats
        assert 'q50' in stats
        assert 'q95' in stats
        assert 'q99' in stats
        assert 'max' in stats

        # Mean should be 3.0 for uniform weights
        assert_allclose(stats['mean'], 3.0, rtol=1e-10)

        # Max should be 5.0
        assert stats['max'] == 5.0


class TestTensorProductMathematicalProperties:
    """Tests verifying mathematical properties of tensor product aggregation."""

    def test_commutativity(self):
        """Order of combination should not matter (approximately)."""
        from survey.aggregation import combine_tensors

        cf1 = np.array([1.0, 2.0, 3.0])
        w1 = np.array([0.2, 0.5, 0.3])
        cf2 = np.array([0.5, 1.5])
        w2 = np.array([0.4, 0.6])

        cf_12, w_12 = combine_tensors(cf1, w1, cf2, w2)
        cf_21, w_21 = combine_tensors(cf2, w2, cf1, w1)

        # Values should be same (possibly different order)
        assert_allclose(sorted(cf_12), sorted(cf_21), rtol=1e-10)

    def test_associativity(self):
        """Grouping of combination should not matter (approximately)."""
        from survey.aggregation import combine_tensors

        cf1 = np.array([1.0, 2.0])
        w1 = np.array([0.5, 0.5])
        cf2 = np.array([0.1, 0.2])
        w2 = np.array([0.5, 0.5])
        cf3 = np.array([0.01])
        w3 = np.array([1.0])

        # (1 + 2) + 3
        cf_12, w_12 = combine_tensors(cf1, w1, cf2, w2)
        cf_123_a, w_123_a = combine_tensors(cf_12, w_12, cf3, w3)

        # 1 + (2 + 3)
        cf_23, w_23 = combine_tensors(cf2, w2, cf3, w3)
        cf_123_b, w_123_b = combine_tensors(cf1, w1, cf_23, w_23)

        # Results should be equivalent
        assert_allclose(sorted(cf_123_a), sorted(cf_123_b), rtol=1e-10)

    def test_mean_additivity(self):
        """Mean of sum should equal sum of means."""
        from survey.aggregation import combine_tensors

        cf1 = np.array([1.0, 2.0, 3.0])
        w1 = np.array([0.2, 0.5, 0.3])
        cf2 = np.array([0.5, 1.0, 1.5])
        w2 = np.array([0.3, 0.4, 0.3])

        mean1 = np.average(cf1, weights=w1)
        mean2 = np.average(cf2, weights=w2)

        cf_combined, w_combined = combine_tensors(cf1, w1, cf2, w2)
        mean_combined = np.average(cf_combined, weights=w_combined)

        assert_allclose(mean_combined, mean1 + mean2, rtol=1e-10)

    def test_variance_additivity(self):
        """Variance of sum of independent RVs equals sum of variances."""
        from survey.aggregation import combine_tensors

        cf1 = np.array([1.0, 2.0, 3.0])
        w1 = np.array([0.2, 0.5, 0.3])
        cf2 = np.array([0.5, 1.0, 1.5])
        w2 = np.array([0.3, 0.4, 0.3])

        mean1 = np.average(cf1, weights=w1)
        var1 = np.average((cf1 - mean1)**2, weights=w1)

        mean2 = np.average(cf2, weights=w2)
        var2 = np.average((cf2 - mean2)**2, weights=w2)

        cf_combined, w_combined = combine_tensors(cf1, w1, cf2, w2)
        mean_combined = np.average(cf_combined, weights=w_combined)
        var_combined = np.average((cf_combined - mean_combined)**2, weights=w_combined)

        # Var(X + Y) = Var(X) + Var(Y) for independent X, Y
        assert_allclose(var_combined, var1 + var2, rtol=1e-10)
