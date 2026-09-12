import numpy as np
import pytest

from mecfs_bio.build_system.task.gwaslab.ldsc_diagnostic import (
    bin_by_ld_score,
    chi2_to_heritability_units,
    fit_line_chi2,
)


def test_bins_are_ordered_by_ld_score_with_equal_counts():
    # Twelve SNPs into three bins -> four SNPs per bin, bins ordered by LD score.
    # Deliberately shuffled input to prove binning is by LD-score rank, not input order.
    ld = np.array([3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0])
    chi2 = np.ones_like(ld)
    bins = bin_by_ld_score(chi2, ld, n_bins=3)
    assert bins.count.tolist() == [4, 4, 4]
    assert bins.mean_ld.tolist() == [1.0, 2.0, 3.0]


def test_bin_means_recover_group_chi2_averages():
    # Sorted so the three equal-count bins are the three contiguous groups; each group's
    # mean chi-square is known by construction.
    ld = np.array([1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0])
    chi2 = np.array([1.0, 3.0, 1.0, 3.0, 2.0, 4.0, 2.0, 4.0, 5.0, 5.0, 7.0, 7.0])
    bins = bin_by_ld_score(chi2, ld, n_bins=3)
    assert bins.mean_chi2.tolist() == [2.0, 3.0, 6.0]


def test_unequal_split_uses_every_snp():
    ld = np.arange(10.0)
    chi2 = np.ones_like(ld)
    bins = bin_by_ld_score(chi2, ld, n_bins=3)
    assert int(bins.count.sum()) == 10
    assert bins.count.max() - bins.count.min() <= 1


def test_se_of_mean_is_sample_std_over_sqrt_count():
    # One bin, two values [1, 5]: sample std (ddof=1) is 2*sqrt(2), so the SE of the mean is
    # 2*sqrt(2)/sqrt(2) = 2.0.
    ld = np.array([1.0, 1.0])
    chi2 = np.array([1.0, 5.0])
    bins = bin_by_ld_score(chi2, ld, n_bins=1)
    assert bins.se_chi2[0] == pytest.approx(2.0)


def test_bin_by_ld_score_rejects_mismatched_lengths():
    with pytest.raises(AssertionError):
        bin_by_ld_score(np.ones(5), np.ones(4), n_bins=2)


def test_heritability_units_rescale_chi2_by_m_over_n():
    # The right-hand axis is chi^2 * M / N, so a heritability slope reads off directly.
    chi2 = np.array([1.0, 2.0, 4.0])
    got = chi2_to_heritability_units(chi2, n=1000.0, m=500.0)
    assert got.tolist() == [0.5, 1.0, 2.0]


def test_fit_line_is_intercept_plus_heritability_slope():
    # E[chi^2] = intercept + (N * h2 / M) * ld ; here slope = 2000 * 0.1 / 1000 = 0.2.
    ld = np.array([0.0, 10.0])
    got = fit_line_chi2(intercept=1.0, h2_obs=0.1, n=2000.0, m=1000.0, ld=ld)
    assert got.tolist() == [1.0, 3.0]
