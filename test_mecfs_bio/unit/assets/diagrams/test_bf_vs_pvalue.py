import numpy as np
from scipy.stats import norm

from mecfs_bio.assets.diagrams.bf_vs_pvalue import curve_points, sampling_variance

_PRIOR_SD_OMEGA = 0.15
_ALLELE_FREQUENCY = 0.3


def test_curve_points_match_density_ratio_and_two_sided_p() -> None:
    """
    The Wakefield ABF is the ratio of the marginal densities of beta_hat under
    the alternative, N(0, V + omega^2), and the null, N(0, V). Checked against
    that ratio and the two-sided p value computed directly, in the range where
    neither underflows.
    """
    beta_hat = 0.05
    sample_size = np.array([1e2, 1e3, 5e3])
    variance = sampling_variance(sample_size, allele_frequency=_ALLELE_FREQUENCY)
    expected_log_bf = norm.logpdf(
        beta_hat, scale=np.sqrt(variance + _PRIOR_SD_OMEGA**2)
    ) - norm.logpdf(beta_hat, scale=np.sqrt(variance))
    expected_p = 2.0 * norm.sf(beta_hat / np.sqrt(variance))

    points = curve_points(
        beta_hat,
        sample_size,
        prior_sd_omega=_PRIOR_SD_OMEGA,
        allele_frequency=_ALLELE_FREQUENCY,
    )

    np.testing.assert_allclose(points.log10_bf, expected_log_bf / np.log(10.0))
    np.testing.assert_allclose(points.neg_log10_p, -np.log10(expected_p))


def test_neg_log10_p_finite_far_past_float_underflow() -> None:
    """
    At z ~ 60 the p value is ~1e-785, below the smallest double; the log survival
    function keeps -log10 p finite.
    """
    points = curve_points(
        0.3,
        np.array([1e5]),
        prior_sd_omega=_PRIOR_SD_OMEGA,
        allele_frequency=_ALLELE_FREQUENCY,
    )
    assert np.all(np.isfinite(points.neg_log10_p))
    assert points.neg_log10_p[0] > 300.0
