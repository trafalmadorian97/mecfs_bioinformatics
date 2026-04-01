import math

import numpy as np
import polars as pl
import pytest

from mecfs_bio.build_system.task.lcv.lcv_core import (
    MomentEstimate,
    as_1d_float_array,
    compute_kappas,
    default_ld_weights,
    estimate_cross_trait_ldsc,
    estimate_mixed_fourth_moments,
    estimate_normalized_scales,
    estimate_trait_ldsc,
    gcp_score_for_value,
    leave_one_block_out_indices,
    validate_equal_length,
    weighted_least_squares,
    weighted_mean,
)

# --- as_1d_float_array ---


def test_as_1d_float_array_converts_list():
    result = as_1d_float_array([1, 2, 3])
    assert result.dtype == np.float64
    np.testing.assert_array_equal(result, [1.0, 2.0, 3.0])


def test_as_1d_float_array_rejects_2d():
    with pytest.raises(ValueError, match="one-dimensional"):
        as_1d_float_array([[1, 2], [3, 4]])


# --- validate_equal_length ---


def test_validate_equal_length_passes():
    validate_equal_length(a=np.array([1, 2]), b=np.array([3, 4]))


def test_validate_equal_length_raises():
    with pytest.raises(ValueError, match="equal length"):
        validate_equal_length(a=np.array([1, 2]), b=np.array([3]))


# --- weighted_mean ---


def test_weighted_mean_uniform_weights():
    assert weighted_mean([1.0, 2.0, 3.0], [1.0, 1.0, 1.0]) == pytest.approx(2.0)


def test_weighted_mean_non_uniform():
    # (1*3 + 2*1) / (3+1) = 5/4
    assert weighted_mean([1.0, 2.0], [3.0, 1.0]) == pytest.approx(5.0 / 4.0)


def test_weighted_mean_single_element():
    assert weighted_mean([7.0], [2.0]) == pytest.approx(7.0)


# --- weighted_least_squares ---


def test_wls_recovers_slope_through_origin():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    y = 3.0 * x
    w = np.ones_like(x)
    beta = weighted_least_squares(x, y, w)
    assert beta[0] == pytest.approx(3.0)


def test_wls_recovers_slope_and_intercept():
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    y = 2.0 * x + 5.0
    design = np.column_stack([x, np.ones_like(x)])
    w = np.ones_like(x)
    beta = weighted_least_squares(design, y, w)
    assert beta[0] == pytest.approx(2.0)
    assert beta[1] == pytest.approx(5.0)


def test_wls_weights_emphasize_subset():
    x = np.array([1.0, 2.0, 3.0])
    y = np.array([10.0, 0.0, 0.0])
    w = np.array([1e8, 1.0, 1.0])
    beta = weighted_least_squares(x, y, w)
    assert beta[0] == pytest.approx(10.0, rel=1e-4)


def test_wls_rejects_3d_design():
    with pytest.raises(ValueError, match="1D or 2D"):
        weighted_least_squares(np.ones((2, 2, 2)), np.ones(2), np.ones(2))


def test_wls_rejects_mismatched_lengths():
    with pytest.raises(ValueError, match="compatible lengths"):
        weighted_least_squares(np.ones(3), np.ones(4), np.ones(3))


# --- default_ld_weights ---


def test_default_ld_weights():
    ld = np.array([0.5, 1.0, 2.0, 10.0])
    w = default_ld_weights(ld)
    # max(1, ld) -> [1, 1, 2, 10], inverse -> [1, 1, 0.5, 0.1]
    np.testing.assert_allclose(w, [1.0, 1.0, 0.5, 0.1])


# --- estimate_trait_ldsc ---


def test_estimate_trait_ldsc_known_slope():
    rng = np.random.default_rng(42)
    n = 5000
    ld = rng.uniform(1, 100, size=n)
    true_slope = 0.5
    true_intercept = 1.0
    chisq = true_slope * ld + true_intercept
    z = np.sqrt(chisq) * rng.choice([-1, 1], size=n)
    w = default_ld_weights(ld)

    result = estimate_trait_ldsc(ld, z, weights=w, chisq_exclude_factor_threshold=30)
    assert result.h2_slope == pytest.approx(true_slope, rel=0.05)
    assert result.intercept == pytest.approx(true_intercept, rel=0.1)


# --- estimate_cross_trait_ldsc ---


def test_estimate_cross_trait_ldsc_runs():
    rng = np.random.default_rng(99)
    n = 3000
    ld = rng.uniform(1, 50, size=n)
    z1 = rng.standard_normal(n) * np.sqrt(0.3 * ld + 1.0)
    z2 = rng.standard_normal(n) * np.sqrt(0.2 * ld + 1.0)
    w = default_ld_weights(ld)

    result = estimate_cross_trait_ldsc(
        ld,
        z1,
        z2,
        weights=w,
        significance_threshold=30,
        h2_slope1=0.3,
        h2_slope2=0.2,
    )
    assert np.isfinite(result.genetic_correlation)
    assert np.isfinite(result.intercept)


# --- estimate_normalized_scales ---


def test_estimate_normalized_scales_positive():
    z1 = np.array([2.0, 3.0, 4.0])
    z2 = np.array([1.5, 2.5, 3.5])
    w = np.ones(3)
    s1, s2 = estimate_normalized_scales(
        z1, z2, weights=w, intercept1=0.0, intercept2=0.0
    )
    assert s1 == pytest.approx(math.sqrt(29.0 / 3.0))
    assert s2 > 0


def test_estimate_normalized_scales_raises_on_negative():
    z1 = np.array([0.5, 0.5, 0.5])
    z2 = np.array([0.5, 0.5, 0.5])
    w = np.ones(3)
    with pytest.raises(ValueError, match="Non-positive"):
        estimate_normalized_scales(z1, z2, weights=w, intercept1=10.0, intercept2=0.0)


# --- estimate_mixed_fourth_moments ---


def test_estimate_mixed_fourth_moments_finite():
    rng = np.random.default_rng(7)
    n = 500
    z1 = rng.standard_normal(n) * 3.0
    z2 = rng.standard_normal(n) * 2.0
    w = np.ones(n)
    k1, k2 = estimate_mixed_fourth_moments(
        z1,
        z2,
        weights=w,
        scale1=3.0,
        scale2=2.0,
        intercept1=1.0,
        intercept2=1.0,
        cross_intercept=0.0,
    )
    assert np.isfinite(k1)
    assert np.isfinite(k2)


# --- leave_one_block_out_indices ---


def test_leave_one_block_out_correct_count():
    blocks = list(leave_one_block_out_indices(100, 5))
    assert len(blocks) == 5


def test_leave_one_block_out_drops_correct_size():
    n_snps, n_blocks = 100, 5
    block_size = n_snps // n_blocks
    for idx_arr in leave_one_block_out_indices(n_snps, n_blocks):
        assert len(idx_arr) == n_snps - block_size


def test_leave_one_block_out_covers_all_indices():
    n_snps, n_blocks = 12, 3
    blocks = list(leave_one_block_out_indices(n_snps, n_blocks))
    counts = np.zeros(n_snps, dtype=int)
    for idx_arr in blocks:
        counts[idx_arr] += 1
    np.testing.assert_array_equal(counts, n_blocks - 1)


def test_leave_one_block_out_rejects_too_few_snps():
    with pytest.raises(ValueError):
        list(leave_one_block_out_indices(3, 10))


def test_leave_one_block_out_rejects_single_block():
    with pytest.raises(ValueError):
        list(leave_one_block_out_indices(10, 1))


# --- compute_kappas ---


def test_compute_kappas_subtracts_3_rho():
    df = pl.DataFrame(
        {
            "rho_g": [0.5, 0.6],
            "mixed_fourth_trait_1": [10.0, 12.0],
            "mixed_fourth_trait_2": [20.0, 22.0],
        }
    )
    rho, k1, k2 = compute_kappas(df)
    np.testing.assert_allclose(rho, [0.5, 0.6])
    np.testing.assert_allclose(k1, [10.0 - 1.5, 12.0 - 1.8])
    np.testing.assert_allclose(k2, [20.0 - 1.5, 22.0 - 1.8])


# --- gcp_score_for_value ---


def test_gcp_score_returns_finite():
    rng = np.random.default_rng(11)
    n = 50
    rho = np.full(n, 0.3)
    kappa1 = rng.standard_normal(n) * 0.1
    kappa2 = rng.standard_normal(n) * 0.1
    stat, discrepancy = gcp_score_for_value(
        0.0,
        rho=rho,
        kappa1=kappa1,
        kappa2=kappa2,
        n_blocks=n,
    )
    assert np.isfinite(stat)
    assert discrepancy.shape == (n,)


def test_gcp_score_near_zero_at_true_gcp():
    """When kappa1/f == f*kappa2, the numerator is zero so the statistic is ~0."""
    rng = np.random.default_rng(0)
    n = 80
    rho = np.full(n, 0.5)
    gcp_true = 0.4
    f = 0.5 ** (-gcp_true)
    c = 1.0
    # Add tiny noise to avoid 0/0 in std
    kappa1 = np.full(n, c * f) + rng.normal(0, 1e-10, n)
    kappa2 = np.full(n, c / f) + rng.normal(0, 1e-10, n)
    stat, _ = gcp_score_for_value(
        gcp_true,
        rho=rho,
        kappa1=kappa1,
        kappa2=kappa2,
        n_blocks=n,
    )
    assert abs(stat) < 1.0


# --- MomentEstimate.as_df ---


def test_moment_estimate_as_df():
    m = MomentEstimate(
        rho_g=0.3,
        mixed_fourth_trait1=1.0,
        mixed_fourth_trait2=2.0,
        cross_intercept=0.05,
        scale1=3.0,
        scale2=4.0,
        intercept1=1.1,
        intercept2=1.2,
    )
    df = m.as_df()
    assert df.shape == (1, 8)
    assert df["rho_g"].item() == pytest.approx(0.3)
    assert df["scale1"].item() == pytest.approx(3.0)
    assert df["mixed_fourth_trait_1"].item() == pytest.approx(1.0)
