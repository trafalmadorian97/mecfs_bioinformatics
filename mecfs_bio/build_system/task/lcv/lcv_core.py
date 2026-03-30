"""
Code implementing the LCV method from
O’Connor, Luke J., and Alkes L. Price. "Distinguishing genetic correlation
from causation across 52 diseases and complex traits."Nature genetics" 50.12 (2018): 1728-1734.

Translated from the original R code using chatgpt, then tweaked
"""
import polars as pl
import pandas as pd
from attrs import frozen
import math
import warnings
from typing import Iterable, Sequence

import numpy as np
import numpy.typing as npt
from scipy.stats import t


FloatArray = npt.NDArray[np.float64]
ArrayLike1D = npt.ArrayLike

@frozen
class TraitLDScoreResult:
    """
    Result of the univariate LDSC-style regression for a single trait.
    """
    h2_slope: float
    intercept: float

@frozen
class CrossTraitLDScoreResult:
    """
    Result of the cross-trait LDSC-style regression.
    """
    genetic_correlation: float
    intercept: float


@frozen
class TraitLDScoreResult:
    """
    Result of the univariate LDSC-style regression for a single trait.
    """
    h2_slope: float
    intercept: float


@frozen
class MomentEstimate:
    """
    Quantities estimated from one dataset (or one jackknife leave-block-out subset).
    """
    rho_g: float
    mixed_fourth_trait1: float
    mixed_fourth_trait2: float
    cross_intercept: float
    scale1: float
    scale2: float
    intercept1: float
    intercept2: float


@frozen
class MomentEstimate:
    """
    Quantities estimated from one dataset (or one jackknife leave-block-out subset).
    """
    rho_g: float
    mixed_fourth_trait1: float
    mixed_fourth_trait2: float
    cross_intercept: float
    scale1: float
    scale2: float
    intercept1: float
    intercept2: float

    def as_df(self)-> pl.DataFrame:
        return pl.DataFrame(
            {
                "rho_g":[self.rho_g],
                "mixed_fourth_trait_1":[self.mixed_fourth_trait1],
                "mixed_fourth_trait_2": [self.mixed_fourth_trait2],
                "cross_intercept":[self.cross_intercept],
                "scale1":[self.scale1],
                "scale2":[self.scale2],
                "intercept1":[self.intercept1],
                "intercept2":[self.intercept2],
            }
        )

@frozen
class JackknifeSummary:
    """
    Block jackknife matrix and the summary quantities derived from it.
    """
    estimates: pl.DataFrame
    rho_est: float
    rho_se: float

@frozen
class LCVResult:
    """
    Final LCV output.
    """
    zscore_gcp_zero: float
    pvalue_gcp_zero_two_sided: float
    posterior_mean_gcp: float
    posterior_se_gcp: float
    rho_est: float
    rho_se: float
    pvalue_gcp_plus_one: float
    pvalue_gcp_minus_one: float
    h2_zscore_trait1: float
    h2_zscore_trait2: float
    gcp_grid: FloatArray
    gcp_weight: FloatArray


def as_1d_float_array(x: ArrayLike1D) -> FloatArray:
    arr = np.asarray(x, dtype=np.float64)
    if arr.ndim != 1:
        raise ValueError(f"array must be one-dimensional")
    return arr

def validate_equal_length(**arrays: FloatArray) -> None:
    lengths = {name: len(arr) for name, arr in arrays.items()}
    if len(set(lengths.values())) != 1:
        msg = ", ".join(f"{k}={v}" for k, v in lengths.items())
        raise ValueError(f"Arrays must have equal length; got {msg}")


def weighted_mean(values: ArrayLike1D, weights: ArrayLike1D) -> float:
    """
    Weighted average:
        sum_i w_i x_i / sum_i w_i
    """
    x = as_1d_float_array(values)
    w = as_1d_float_array(weights)
    validate_equal_length(values=x, weights=w)
    return float(np.sum(w * x) / np.sum(w))

def weighted_least_squares(design: npt.ArrayLike, response: ArrayLike1D, weights: ArrayLike1D) -> FloatArray:
    """
    Solve the weighted least squares normal equations:

        beta = (X^T W X)^(-1) X^T W y

    Parameters
    ----------
    design:
        Shape (n,) or (n, p). A 1D input is treated as a single predictor column.
    response:
        Shape (n,)
    weights:
        Shape (n,)
    """
    y = as_1d_float_array(response)
    w = as_1d_float_array(weights)

    X = np.asarray(design, dtype=np.float64)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    elif X.ndim != 2:
        raise ValueError("design must be 1D or 2D")

    if X.shape[0] != len(y) or len(y) != len(w):
        raise ValueError("design, response, and weights must have compatible lengths")

    WX = X * w[:, None]
    lhs = X.T @ WX
    rhs = X.T @ (w * y)
    return np.linalg.solve(lhs, rhs)


def default_ld_weights(ld_scores: ArrayLike1D) -> FloatArray:
    ell = as_1d_float_array(ld_scores)
    return 1.0 / np.maximum(1.0, ell)


def estimate_trait_ldsc(
    ld_scores: ArrayLike1D,
    z_scores: ArrayLike1D,
    *,
    weights: ArrayLike1D,
    significance_threshold: float,
) -> TraitLDScoreResult:
    """
    Estimate the trait-specific LDSC-style slope/intercept from z^2 ~ ld.

    """
    ld = as_1d_float_array(ld_scores)
    z = as_1d_float_array(z_scores)
    w = as_1d_float_array(weights)
    validate_equal_length(ld_scores=ld, z_scores=z, weights=w)


    chisq = z**2
    mean_chisq = float(np.mean(chisq))
    keep = chisq <= significance_threshold * mean_chisq
    chisq_keep = chisq[keep]

    ld_score_and_constant_keep = np.concat([ld[keep].reshape((-1,1)), np.ones(np.sum(keep), dtype=np.float64).reshape((-1,1))], axis=1)
    regression_results = weighted_least_squares(ld_score_and_constant_keep, chisq_keep, w[keep])
    intercept = float(regression_results[1])

    slope = float(weighted_least_squares(ld, chisq - intercept, w)[0])

    return TraitLDScoreResult(h2_slope=slope, intercept=intercept)


def estimate_cross_trait_ldsc(
    ld_scores: ArrayLike1D,
    z1: ArrayLike1D,
    z2: ArrayLike1D,
    *,
    weights: ArrayLike1D,
    significance_threshold: float,
    h2_slope1: float,
    h2_slope2: float,
) -> CrossTraitLDScoreResult:
    """
    Estimate the cross-trait LDSC-style slope from z1*z2 ~ ld_arr, then normalize
    by sqrt(h2_slope1 * h2_slope2) to obtain the genetic correlation estimate.
    """
    ld_arr = as_1d_float_array(ld_scores)
    z1_arr = as_1d_float_array(z1)
    z2_arr = as_1d_float_array(z2)
    w_arr = as_1d_float_array(weights)
    validate_equal_length(ld_scores=ld_arr, z1=z1_arr, z2=z2_arr, weights=w_arr)


    mean_chisq1 = float(np.mean(z1_arr**2))
    mean_chisq2 = float(np.mean(z2_arr**2))
    keep = (
        (z1_arr**2 < significance_threshold * mean_chisq1)
        & (z2_arr**2 < significance_threshold * mean_chisq2)
    )

    ld_score_and_constant_keep = np.concat([ld_arr[keep].reshape(-1,1), np.ones(np.sum(keep).reshape(-1,1), dtype=np.float64)], axis=1)
    regression_results = weighted_least_squares(ld_score_and_constant_keep, z1_arr[keep] * z2_arr[keep], w_arr[keep])
    cross_intercept = float(regression_results[1])

    slope = float(weighted_least_squares(ld_arr, z1_arr * z2_arr - cross_intercept, w_arr)[0])
    rho = slope / math.sqrt(h2_slope1 * h2_slope2)
    return CrossTraitLDScoreResult(genetic_correlation=rho, intercept=cross_intercept)

def estimate_normalized_scales(
    z1: ArrayLike1D,
    z2: ArrayLike1D,
    *,
    weights: ArrayLike1D,
    intercept1: float,
    intercept2: float,
) -> tuple[float, float]:
    """
    Compute the normalization scales used by LCV:
        s1 = sqrt(E_w[z1^2] - intercept1)
        s2 = sqrt(E_w[z2^2] - intercept2)

    These values are used to normalize z scores
    """
    z1_arr = as_1d_float_array(z1)
    z2_arr = as_1d_float_array(z2)
    w = as_1d_float_array(weights)
    validate_equal_length(z1=z1_arr, z2=z2_arr, weights=w)

    s1_sq = weighted_mean(z1_arr**2, w) - intercept1
    s2_sq = weighted_mean(z2_arr**2, w) - intercept2

    if s1_sq <= 0 or s2_sq <= 0:
        raise ValueError(
            "Non-positive normalization scale encountered. "
            "This often reflects unstable or negative heritability estimates."
        )

    return math.sqrt(s1_sq), math.sqrt(s2_sq)


def estimate_mixed_fourth_moments(
    z1: ArrayLike1D,
    z2: ArrayLike1D,
    *,
    weights: ArrayLike1D,
    scale1: float,
    scale2: float,
    intercept1: float,
    intercept2: float,
    cross_intercept: float,
) -> tuple[float, float]:
    """
    Estimate the two corrected mixed fourth moments used by LCV.

    These correspond to the R variables:
        k41, k42

    My note: the computations of kappa1 and kappa2 below come from equation 6 in the original paper.

    """
    z1_arr = as_1d_float_array(z1)
    z2_arr = as_1d_float_array(z2)
    w = as_1d_float_array(weights)
    validate_equal_length(z1=z1_arr, z2=z2_arr, weights=w)

    nz1 = z1_arr / scale1
    nz2 = z2_arr / scale2

    kappa1_expression = (
        nz2 * nz1**3
        - 3.0 * nz1 * nz2 * (intercept1 / scale1**2)
        - 3.0 * (nz1**2 - intercept1 / scale1**2) * cross_intercept / (scale1 * scale2)
    )
    kappa2_expression = (
        nz1 * nz2**3
        - 3.0 * nz1 * nz2 * (intercept2 / scale2**2)
        - 3.0 * (nz2**2 - intercept2 / scale2**2) * cross_intercept / (scale1 * scale2)
    )

    kappa1 = weighted_mean(kappa1_expression, w)
    kappa2 = weighted_mean(kappa2_expression, w)
    return kappa1, kappa2



def estimate_lcv_moments(
    ld_scores: ArrayLike1D,
    z1: ArrayLike1D,
    z2: ArrayLike1D,
    *,
    significance_threshold: float ,
    weights: ArrayLike1D | None = None,
) -> MomentEstimate:
    """
    One-shot estimation of all LCV moments on a given SNP set.

    This is the refactored counterpart of EstimateK4(..., nargout=8).
    """
    ld_arr = as_1d_float_array(ld_scores)
    z1_arr = as_1d_float_array(z1)
    z2_arr = as_1d_float_array(z2)
    validate_equal_length(ld_scores=ld_arr, z1=z1_arr, z2=z2_arr)

    w = default_ld_weights(ld_arr) if weights is None else as_1d_float_array(weights)
    validate_equal_length(ld_scores=ld_arr, weights=w)

    trait1 = estimate_trait_ldsc(
        ld_arr,
        z1_arr,
        weights=w,
        significance_threshold=significance_threshold,
    )
    trait2 = estimate_trait_ldsc(
        ld_arr,
        z2_arr,
        weights=w,
        significance_threshold=significance_threshold,
    )

    cross = estimate_cross_trait_ldsc(
        ld_arr,
        z1_arr,
        z2_arr,
        weights=w,
        significance_threshold=significance_threshold,
        h2_slope1=trait1.h2_slope,
        h2_slope2=trait2.h2_slope,
    )

    scale1, scale2 = estimate_normalized_scales(
        z1_arr,
        z2_arr,
        weights=w,
        intercept1=trait1.intercept,
        intercept2=trait2.intercept,
    )

    kappa1, kappa2 = estimate_mixed_fourth_moments(
        z1_arr,
        z2_arr,
        weights=w,
        scale1=scale1,
        scale2=scale2,
        intercept1=trait1.intercept,
        intercept2=trait2.intercept,
        cross_intercept=cross.intercept,
    )

    return MomentEstimate(
        rho_g=cross.genetic_correlation,
        mixed_fourth_trait1=kappa1,
        mixed_fourth_trait2=kappa2,
        cross_intercept=cross.intercept,
        scale1=scale1,
        scale2=scale2,
        intercept1=trait1.intercept,
        intercept2=trait2.intercept,
    )


def leave_one_block_out_indices(n_snps: int, n_blocks: int) -> Iterable[FloatArray]:
    """
    Yield index arrays corresponding to leave-one-block-out subsets.

    Blocks are consecutive in the current SNP ordering, matching the original R code.
    """
    if n_blocks < 2:
        raise ValueError("n_blocks must be at least 2")
    if n_snps < n_blocks:
        raise ValueError("n_snps must be at least n_blocks")

    block_size = n_snps // n_blocks
    if block_size == 0:
        raise ValueError("block size is zero; reduce n_blocks")

    for block in range(n_blocks):
        if block == 0:
            yield np.arange(block_size, n_snps)
        elif block == n_blocks - 1:
            yield np.arange(0, block_size * block)
        else:
            left = np.arange(0, block * block_size)
            right = np.arange((block + 1) * block_size, n_snps)
            yield np.concatenate([left, right])


def compute_jackknife_summary(
    ld_scores: ArrayLike1D,
    z1: ArrayLike1D,
    z2: ArrayLike1D,
    *,
    significance_threshold: float ,
    n_blocks: int = 100,
    weights: ArrayLike1D | None = None,
) -> JackknifeSummary:
    """
    Compute leave-one-block-out estimates of the eight LCV moment quantities.
    """
    ld_arr = as_1d_float_array(ld_scores)
    z1_arr = as_1d_float_array(z1)
    z2_arr = as_1d_float_array(z2)
    validate_equal_length(ld_scores=ld_arr, z1=z1_arr, z2=z2_arr)

    w = default_ld_weights(ld_arr) if weights is None else as_1d_float_array(weights)
    validate_equal_length(ld_scores=ld_arr, weights=w)

    estimates = []

    for i, keep_idx in enumerate(leave_one_block_out_indices(len(ld_arr), n_blocks)):
        estimate = estimate_lcv_moments(
            ld_arr[keep_idx],
            z1_arr[keep_idx],
            z2_arr[keep_idx],
            weights=w[keep_idx],
            significance_threshold=significance_threshold,
        )
        estimates.append(estimate.as_df())

    if np.isnan(estimates).any():
        raise ValueError(
            "NaNs produced during jackknife estimation. "
            "This often indicates negative/unstable heritability estimates or misordered SNPs."
        )
    estimates_df =pl.concat(estimates,how="vertical").collect()

    rho_est = float(np.mean(estimates_df["rho_g"].to_numpy()))
    rho_se = float(np.std(estimates_df["rho_g"].to_numpy(), ddof=1) * math.sqrt(n_blocks + 1))

    return JackknifeSummary(estimates=estimates_df, rho_est=rho_est, rho_se=rho_se)


def compute_kappas(jackknife_estimates: pl.DataFrame) -> tuple[FloatArray, FloatArray, FloatArray]:
    """
    Extract rho, kappa1, kappa2 and subtract the Gaussian null term 3*rho
    from the fourth moments, matching the original implementation.

    My notes:
    From the equation at the top of page 11
    E(\alpha_1^3\alpha_2)=q_1^3q_2*(kurtosis_1) + 3q_1q_2
    =q_1^3q_2*(kurtosis_pi) + 3*rho_g

    The key step of this function is to subtract 3*rho_g from estimates of E(\alpha_1^3\alpha_2)  and so return
    the kappas, which are estimates of q_1^3q_2*(kurtosis_pi)
    """
    rho = jackknife_estimates["rho_g"].to_numpy()
    kappa1 = jackknife_estimates["mixed_fourth_trait_1"].to_numpy() - 3.0 * rho
    kappa2 = jackknife_estimates["mixed_fourth_trait_2"].to_numpy() - 3.0 * rho
    return rho, kappa1, kappa2


def gcp_score_for_value(
    gcp: float,
    *,
    rho: FloatArray,
    kappa1: FloatArray,
    kappa2: FloatArray,
    n_blocks: int,
) -> tuple[float, FloatArray]:
    """
    Compute the jackknife-based discrepancy statistic for a candidate gcp value.

    This is a direct but renamed version of the inner loop in RunLCV.R.

    My note:Need to investigate this.  I think the idea is that for the true gcp, numerator=0.  But how do we show this?

    """
    scale_factor = np.abs(rho) ** (-gcp)

    numerator = kappa1 / scale_factor - scale_factor * kappa2
    denominator = np.maximum(
        1.0 / np.abs(rho),
        np.sqrt(kappa1**2 / scale_factor**2 + kappa2**2 * scale_factor**2),
    )
    standardized_discrepancy = numerator / denominator

    statistic = (
        np.mean(standardized_discrepancy)
        / np.std(standardized_discrepancy, ddof=1)
        / math.sqrt(n_blocks + 1)
    )
    return float(statistic), standardized_discrepancy