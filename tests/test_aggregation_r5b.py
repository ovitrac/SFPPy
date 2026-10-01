"""
Aggregation primitives (survey.aggregation).

Fast unit tests on synthetic data (no Patankar, no PubChem):

  sources_share_t_axis      : invariant check, raises on mismatch
  sum_sources_at_fixed_t    : Cp0 tensor product at a fixed t cell
  combine_sources_shared_t  : full group aggregation with shared t
  jaccard_similarity        : diagnostic flag

Integration cases:
  - single source (pass-through)
  - two sources, 1-step, both active
  - two sources, 2-step, both active (shared t1 AND t2)
  - two sources, mixed step-2 activity (one active, one removed at
    step 2) — inactive source contributes constant CF across t2
  - linearity: Cp0 tensor product preserves superposition

@project: SFPPy — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from survey.aggregation import (  # noqa: E402
    sources_share_t_axis,
    sum_sources_at_fixed_t,
    combine_sources_shared_t,
    jaccard_similarity,
)


# ======================================================================
# Synthetic source factories
# ======================================================================

def _src_1step(cf_tensor, time_vals, time_weights, conc_vals, conc_weights):
    return {
        "CF_tensor": np.asarray(cf_tensor, dtype=float),
        "time_vals": np.asarray(time_vals, dtype=float),
        "time_weights": np.asarray(time_weights, dtype=float),
        "conc_vals": np.asarray(conc_vals, dtype=float),
        "conc_weights": np.asarray(conc_weights, dtype=float),
        # step2_active absent → 1-step scenario
    }


def _src_2step(cf_tensor, time_vals, time_weights,
               time2_vals, time2_weights, conc_vals, conc_weights):
    return {
        "CF_tensor": np.asarray(cf_tensor, dtype=float),
        "time_vals": np.asarray(time_vals, dtype=float),
        "time_weights": np.asarray(time_weights, dtype=float),
        "time2_vals": np.asarray(time2_vals, dtype=float),
        "time2_weights": np.asarray(time2_weights, dtype=float),
        "conc_vals": np.asarray(conc_vals, dtype=float),
        "conc_weights": np.asarray(conc_weights, dtype=float),
        "step2_active": True,
    }


# ======================================================================
# sources_share_t_axis
# ======================================================================

def test_share_t_axis_empty_raises():
    with pytest.raises(ValueError, match="non-empty"):
        sources_share_t_axis([])


def test_share_t_axis_ok_1step():
    t1 = np.array([1.0, 2.0, 3.0])
    w = np.array([0.2, 0.5, 0.3])
    s1 = _src_1step(np.zeros((3, 2)), t1, w, [1., 2.], [0.5, 0.5])
    s2 = _src_1step(np.zeros((3, 4)), t1, w, [1., 2., 3., 4.], [0.25]*4)
    assert sources_share_t_axis([s1, s2]) is True


def test_share_t_axis_mismatch_raises():
    t1_a = np.array([1.0, 2.0])
    t1_b = np.array([1.0, 3.0])
    w = np.array([0.5, 0.5])
    s1 = _src_1step(np.zeros((2, 2)), t1_a, w, [1., 2.], [0.5, 0.5])
    s2 = _src_1step(np.zeros((2, 2)), t1_b, w, [1., 2.], [0.5, 0.5])
    with pytest.raises(ValueError, match="time_vals mismatch"):
        sources_share_t_axis([s1, s2])


def test_share_t_axis_ok_2step():
    t1 = np.array([1.0, 2.0])
    w1 = np.array([0.6, 0.4])
    t2 = np.array([10.0, 20.0, 30.0])
    w2 = np.array([0.2, 0.5, 0.3])
    s1 = _src_2step(np.zeros((2, 3, 2)), t1, w1, t2, w2, [1., 2.], [0.5, 0.5])
    s2 = _src_2step(np.zeros((2, 3, 4)), t1, w1, t2, w2,
                     [1., 2., 3., 4.], [0.25]*4)
    assert sources_share_t_axis([s1, s2]) is True


def test_share_t_axis_mixed_step2_only_active_share_t2():
    t1 = np.array([1.0, 2.0])
    w1 = np.array([0.6, 0.4])
    t2 = np.array([10.0, 20.0, 30.0])
    w2 = np.array([0.2, 0.5, 0.3])
    active = _src_2step(np.zeros((2, 3, 2)), t1, w1, t2, w2,
                         [1., 2.], [0.5, 0.5])
    # Inactive source has a 2-D CF_tensor and no time2_*; it shares t1
    # only. Must NOT raise.
    inactive = _src_1step(np.zeros((2, 2)), t1, w1, [1., 2.], [0.5, 0.5])
    assert sources_share_t_axis([active, inactive]) is True


# ======================================================================
# sum_sources_at_fixed_t — Cp0 tensor product at fixed t
# ======================================================================

def test_sum_at_fixed_t_single_source_passthrough():
    """One source → the Cp0 tensor product is a no-op."""
    cf = np.array([[1.0, 2.0, 3.0],
                    [10.0, 20.0, 30.0]])  # (n_t, n_cp0)
    s = _src_1step(cf, [1., 2.], [0.5, 0.5],
                     [0.1, 1.0, 10.0], [0.2, 0.5, 0.3])
    cf_out, w_out = sum_sources_at_fixed_t([s], t1_idx=1)
    assert np.allclose(cf_out, [10.0, 20.0, 30.0])
    assert np.allclose(w_out, [0.2, 0.5, 0.3])


def test_sum_at_fixed_t_two_sources_cp0_tensor_product():
    """Two sources → n_cp0_1 * n_cp0_2 cells via tensor product."""
    cf1 = np.array([[10.0, 20.0]])          # shape (1, 2)
    cf2 = np.array([[100.0, 200.0, 300.0]]) # shape (1, 3)
    s1 = _src_1step(cf1, [1.0], [1.0], [0.1, 0.2], [0.3, 0.7])
    s2 = _src_1step(cf2, [1.0], [1.0], [1, 2, 3], [0.1, 0.6, 0.3])
    cf_out, w_out = sum_sources_at_fixed_t([s1, s2], t1_idx=0)
    # Expected: CF_1 + CF_2, outer-weight = w_1 * w_2
    # cf_sums[i,j] = CF_1[i] + CF_2[j]
    expected_cf = np.array([
        [10+100, 10+200, 10+300],
        [20+100, 20+200, 20+300],
    ]).ravel()
    expected_w = np.outer([0.3, 0.7], [0.1, 0.6, 0.3]).ravel()
    # Order-agnostic comparison
    order_out = np.argsort(cf_out)
    order_exp = np.argsort(expected_cf)
    assert np.allclose(cf_out[order_out], expected_cf[order_exp])
    assert np.allclose(w_out[order_out], expected_w[order_exp])


def test_sum_at_fixed_t_step2_inactive_slices_t1_only():
    """Inactive source's 2-D CF_tensor is sliced by t1_idx only."""
    cf_active = np.array([[[1.0, 2.0], [10.0, 20.0], [100.0, 200.0]]])  # (1, 3, 2)
    cf_inactive = np.array([[5.0, 50.0]])                                # (1, 2)
    active = _src_2step(cf_active, [1.0], [1.0], [1., 2., 3.], [0.25, 0.5, 0.25],
                         [0.1, 1.0], [0.5, 0.5])
    inactive = _src_1step(cf_inactive, [1.0], [1.0], [0.1, 1.0], [0.5, 0.5])
    # At t2_idx=1 (mid): active gives [10, 20]; inactive gives [5, 50].
    cf_out, w_out = sum_sources_at_fixed_t(
        [active, inactive], t1_idx=0, t2_idx=1)
    expected = np.array([10+5, 10+50, 20+5, 20+50])
    order_out = np.argsort(cf_out)
    order_exp = np.argsort(expected)
    assert np.allclose(cf_out[order_out], expected[order_exp])


# ======================================================================
# combine_sources_shared_t — full pipeline
# ======================================================================

def test_combine_single_source_bit_identical_to_marginal():
    """Single-source group: output equals the source's marginal CF × Cp0."""
    cf = np.array([[10.0, 20.0, 30.0],
                    [100.0, 200.0, 300.0]])
    t_w = np.array([0.3, 0.7])
    c_w = np.array([0.2, 0.5, 0.3])
    s = _src_1step(cf, [1., 2.], t_w, [0.1, 1., 10.], c_w)
    cf_out, w_out = combine_sources_shared_t([s])
    expected_w = np.outer(t_w, c_w).ravel()
    # Normalised
    assert abs(w_out.sum() - 1.0) < 1e-12
    assert abs(expected_w.sum() - 1.0) < 1e-12

    # Compare means (moment check)
    m_out = np.average(cf_out, weights=w_out)
    m_expected = np.average(cf.ravel(), weights=expected_w)
    assert abs(m_out - m_expected) < 1e-12 * max(abs(m_expected), 1.0)


def test_combine_two_sources_1step_shared_t():
    """Two 1-step sources, shared t1: Cp0 tensor product + shared t1 sum."""
    # Source 1: CF = t (constant across Cp0 for simplicity)
    cf1 = np.array([[1.0], [2.0], [3.0]])  # (n_t=3, n_cp0=1)
    # Source 2: CF = 10*t
    cf2 = np.array([[10.0], [20.0], [30.0]])
    t = np.array([1., 2., 3.])
    tw = np.array([0.25, 0.5, 0.25])
    s1 = _src_1step(cf1, t, tw, [1.0], [1.0])
    s2 = _src_1step(cf2, t, tw, [1.0], [1.0])

    cf_out, w_out = combine_sources_shared_t([s1, s2])
    # Since each source has n_cp0=1 (degenerate), the Cp0 tensor has one cell.
    # Sum at each t: CF_1(t) + CF_2(t) = 11*t. Shared t weight: tw[i].
    expected_cf = np.array([11.0, 22.0, 33.0])
    expected_w = tw
    # Same order
    order_out = np.argsort(cf_out)
    order_exp = np.argsort(expected_cf)
    assert np.allclose(cf_out[order_out], expected_cf[order_exp])
    assert np.allclose(w_out[order_out], expected_w[order_exp])


def test_combine_2step_both_active():
    """Two 2-step sources, shared (t1, t2): 3-D sum over all cells."""
    # Source 1: CF = t1 + 0*t2 + 0*Cp0 (just t1)
    cf1 = np.array([[[1.0], [1.0]], [[2.0], [2.0]]])   # (n_t1=2, n_t2=2, n_cp0=1)
    cf2 = np.array([[[10.0], [10.0]], [[20.0], [20.0]]])
    t1 = np.array([1., 2.])
    t1w = np.array([0.6, 0.4])
    t2 = np.array([10., 20.])
    t2w = np.array([0.3, 0.7])
    s1 = _src_2step(cf1, t1, t1w, t2, t2w, [1.0], [1.0])
    s2 = _src_2step(cf2, t1, t1w, t2, t2w, [1.0], [1.0])

    cf_out, w_out = combine_sources_shared_t([s1, s2])
    # At each (t1_i, t2_j): CF_sum = CF_1[i,j] + CF_2[i,j] = 11*t1_i.
    # Marginalised weights: outer(t1w, t2w).ravel().
    # Cells: n_t1 * n_t2 = 4.
    assert len(cf_out) == 4
    order = np.argsort(cf_out)
    expected_cf = np.array([11.0, 11.0, 22.0, 22.0])
    expected_w = np.array([
        t1w[0]*t2w[0], t1w[0]*t2w[1], t1w[1]*t2w[0], t1w[1]*t2w[1]
    ])
    # sort expected by cf
    exp_order = np.argsort(expected_cf)
    assert np.allclose(cf_out[order], expected_cf[exp_order])
    assert np.allclose(w_out[order], expected_w[exp_order])


def test_combine_mixed_step2_activity():
    """
    One step-2-active source + one step-2-inactive source.
    Inactive source contributes constant CF regardless of t2.
    """
    # Active source: CF(t1, t2) = t2 (ignores t1)
    cf_a = np.array([[[1.0], [5.0], [10.0]],
                      [[1.0], [5.0], [10.0]]])   # (n_t1=2, n_t2=3, n_cp0=1)
    # Inactive source: CF(t1) = 100 (constant across t1 for simplicity)
    cf_i = np.array([[100.0], [100.0]])          # (n_t1=2, n_cp0=1)
    t1 = np.array([1., 2.])
    t1w = np.array([0.5, 0.5])
    t2 = np.array([10., 20., 30.])
    t2w = np.array([0.25, 0.5, 0.25])
    active = _src_2step(cf_a, t1, t1w, t2, t2w, [1.0], [1.0])
    inactive = _src_1step(cf_i, t1, t1w, [1.0], [1.0])

    cf_out, w_out = combine_sources_shared_t([active, inactive])
    # At each (t1_i, t2_j): CF = CF_active[i,j] + 100 = [101, 105, 110]
    # (same for each t1). n_cells = 2 * 3 = 6.
    assert len(cf_out) == 6
    # weights normalised to 1
    assert abs(w_out.sum() - 1.0) < 1e-12
    # All cells should have CF in {101, 105, 110}
    uniq = np.unique(np.round(cf_out, 6))
    assert np.allclose(sorted(uniq), [101.0, 105.0, 110.0])


def test_combine_weights_sum_to_one():
    """Normalisation invariant: Σ w == 1 for any non-empty group."""
    t = np.array([1., 2., 3.])
    tw = np.array([0.2, 0.5, 0.3])
    c = np.array([1., 10.])
    cw = np.array([0.3, 0.7])
    cf = np.random.default_rng(42).uniform(0, 10, size=(3, 2))
    s1 = _src_1step(cf, t, tw, c, cw)
    s2 = _src_1step(cf * 2, t, tw, c, cw)
    _, w = combine_sources_shared_t([s1, s2])
    assert abs(w.sum() - 1.0) < 1e-12


def test_combine_linearity_in_cp0():
    """
    Scaling one source's Cp0 grid by α scales the post-aggregation CF
    distribution's mean by α (in the limit of that source dominating).
    A weak but useful invariant check.
    """
    t = np.array([1.0])
    tw = np.array([1.0])
    c = np.array([1., 10., 100.])
    cw = np.array([0.1, 0.6, 0.3])
    # Source 1: CF = Cp0 (identity)
    cf1 = np.array([c])                    # (1, 3)
    # Source 2: CF = 0 (contributes 0)
    cf2 = np.zeros((1, 3))
    s1 = _src_1step(cf1, t, tw, c, cw)
    s2 = _src_1step(cf2, t, tw, c, cw)
    cf_out, w_out = combine_sources_shared_t([s1, s2])
    mean = np.average(cf_out, weights=w_out)
    expected = np.average(c, weights=cw)   # E[Cp0]
    assert abs(mean - expected) < 1e-12 * abs(expected)


# ======================================================================
# jaccard_similarity
# ======================================================================

def test_jaccard_singleton_is_one():
    assert jaccard_similarity([{"A", "B"}]) == 1.0


def test_jaccard_identical_sets():
    assert jaccard_similarity([{"A", "B"}, {"A", "B"}]) == 1.0


def test_jaccard_disjoint():
    assert jaccard_similarity([{"A"}, {"B"}]) == 0.0


def test_jaccard_partial_overlap():
    j = jaccard_similarity([{"A", "B", "C"}, {"B", "C", "D"}])
    assert abs(j - 2/4) < 1e-12


def test_jaccard_min_over_all_pairs():
    # Three sets: AB, AC, BC → pairwise 1/3, 1/3, 1/3.
    j = jaccard_similarity([{"A", "B"}, {"A", "C"}, {"B", "C"}])
    assert abs(j - 1/3) < 1e-12


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
