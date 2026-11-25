"""Weighted ridge regression accumulated per chromosome block.

Pure numpy: no polars, parquet, or Task imports. A ChromRidgeBlock holds the
weighted cross-product sufficient statistics of one chromosome; combining blocks
lets any leave-one-chromosome-out (or odd/even) fit be formed as a sum, so the
full design matrix is never held in memory. Ridge is fit on centered, per-column
standardized features and unstandardized back; weights of 1 recover ordinary
ridge.

With features x_i (length p), target y_i and weights w_i, the fitted model is
y_i ~ intercept + x_i . beta_raw, where beta minimizes

    sum_i w_i (y_i - mean_y - z_i . beta_std)^2 + alpha ||beta_std||^2

over the standardized features z_ij = (x_ij - mean_j) / sd_j.
"""

from collections.abc import Mapping, Sequence

import numpy as np
from attrs import frozen


@frozen(slots=True)
class ChromRidgeBlock:
    """Weighted cross-product sufficient statistics for one chromosome.

    For the SNPs i in one chromosome, let x_i be the length-p feature row, y_i the
    target, and w_i >= 0 the regression weight. Every field is a weighted sum over
    those SNPs:

        sw   = sum_i w_i                  scalar (total weight; = n when w == 1)
        swx  = sum_i w_i x_i              length p (weighted feature sum)
        swxx = sum_i w_i x_i x_i^T        p x p (weighted Gram, X^T W X)
        swxy = sum_i w_i x_i y_i          length p (weighted feature-target, X^T W y)
        swy  = sum_i w_i y_i              scalar
        swyy = sum_i w_i y_i^2            scalar

    These are exactly the statistics the weighted-ridge normal equations need, and
    every one is additive across chromosomes, so any leave-one-out or odd/even
    training set is formed by summing the relevant blocks (see combine).
    """

    sw: float
    swx: np.ndarray
    swxx: np.ndarray
    swxy: np.ndarray
    swy: float
    swyy: float

    def __attrs_post_init__(self) -> None:
        p = self.swx.shape[0]
        assert self.swx.shape == (p,), f"swx must be 1-d, got {self.swx.shape}"
        assert self.swxx.shape == (p, p), (
            f"swxx must be ({p}, {p}), got {self.swxx.shape}"
        )
        assert self.swxy.shape == (p,), f"swxy must be ({p},), got {self.swxy.shape}"
        for arr in (self.swx, self.swxx, self.swxy):
            assert arr.dtype.kind == "f", f"expected float array, got {arr.dtype}"
        assert self.sw > 0, f"block total weight must be positive, got {self.sw}"


@frozen(slots=True)
class StandardizedSystem:
    """The weighted ridge system after centering and per-column standardization.

    Built from a ChromRidgeBlock. With the standardized feature
    z_ij = (x_ij - mean_j) / sd_j, the fields are:

        mean    length p  weighted per-feature mean:
                          mean_j = (sum_i w_i x_ij) / (sum_i w_i) = swx_j / sw
        sd      length p  weighted per-feature standard deviation (zeros replaced by 1):
                          sd_j = sqrt(swxx_jj / sw - mean_j^2)
        mean_y  scalar    weighted mean of the target: mean_y = swy / sw
        g_std   p x p     weighted Gram of the standardized features:
                          g_std[j,k] = sum_i w_i z_ij z_ik
                                     = (swxx - sw * outer(mean, mean))[j,k] / (sd_j sd_k)
        b_std   length p  weighted standardized feature-target cross term:
                          b_std[j] = sum_i w_i z_ij (y_i - mean_y)
                                   = (swxy - mean * swy)[j] / sd_j
                          (the y-centering term vanishes because sum_i w_i z_ij = 0)

    solve returns beta_std from (g_std + alpha I) beta_std = b_std; fit then
    unstandardizes: beta_raw = beta_std / sd, intercept = mean_y - beta_raw . mean.
    """

    g_std: np.ndarray
    b_std: np.ndarray
    mean: np.ndarray
    sd: np.ndarray
    mean_y: float

    def __attrs_post_init__(self) -> None:
        p = self.b_std.shape[0]
        assert self.g_std.shape == (p, p), f"g_std must be ({p}, {p})"
        assert self.b_std.shape == (p,)
        assert self.mean.shape == (p,)
        assert self.sd.shape == (p,)
        assert np.all(self.sd > 0), "sd must be strictly positive"


@frozen(slots=True)
class RidgeFit:
    """A fitted ridge model: y ~ intercept + x . beta_raw.

    beta_std is the same coefficient vector on the standardized feature scale
    (beta_std = beta_raw * sd), comparable across features.
    """

    beta_raw: np.ndarray
    beta_std: np.ndarray
    intercept: float

    def __attrs_post_init__(self) -> None:
        assert self.beta_raw.ndim == 1
        assert self.beta_std.shape == self.beta_raw.shape


@frozen(slots=True)
class AlphaSelection:
    """The ridge penalty chosen by leave-one-chromosome-out cross-validation,
    with its mean held-out weighted R^2, the per-held-out-chromosome R^2, and the
    mean held-out R^2 of every candidate penalty (the cross-validation curve)."""

    alpha: float
    mean_r2: float
    r2_per_chrom: dict[int, float]
    mean_r2_by_alpha: dict[float, float]

    def __attrs_post_init__(self) -> None:
        assert self.mean_r2_by_alpha[self.alpha] == self.mean_r2, (
            "the selected alpha's score must match its entry in the curve"
        )


def accumulate_block(
    x: np.ndarray, y: np.ndarray, w: np.ndarray | None = None
) -> ChromRidgeBlock:
    """Compute one block's weighted sufficient statistics from an (n, p) feature
    matrix x, a length-n target y, and optional length-n weights w (default 1)."""
    assert x.ndim == 2, f"x must be 2-d, got shape {x.shape}"
    assert y.shape == (x.shape[0],), f"y must have shape ({x.shape[0]},), got {y.shape}"
    x = x.astype(np.float64)
    y = y.astype(np.float64)
    if w is None:
        w = np.ones(x.shape[0])
    assert w.shape == y.shape, f"w must have shape {y.shape}, got {w.shape}"
    assert np.all(w >= 0), "weights must be non-negative"
    w = w.astype(np.float64)
    wx = x * w[:, None]
    return ChromRidgeBlock(
        sw=float(w.sum()),
        swx=wx.sum(axis=0),
        swxx=x.T @ wx,
        swxy=wx.T @ y,
        swy=float((w * y).sum()),
        swyy=float((w * y * y).sum()),
    )


def combine(blocks: Sequence[ChromRidgeBlock]) -> ChromRidgeBlock:
    """Sum blocks into the block of their union (all fields are additive)."""
    assert len(blocks) > 0, "cannot combine zero blocks"
    return ChromRidgeBlock(
        sw=sum(b.sw for b in blocks),
        swx=np.sum([b.swx for b in blocks], axis=0),
        swxx=np.sum([b.swxx for b in blocks], axis=0),
        swxy=np.sum([b.swxy for b in blocks], axis=0),
        swy=sum(b.swy for b in blocks),
        swyy=sum(b.swyy for b in blocks),
    )


def standardized_system(block: ChromRidgeBlock) -> StandardizedSystem:
    """Center and standardize a block's normal equations (see StandardizedSystem).

    The centered Gram follows from expanding the outer product:
    sum_i w_i (x_i - mean)(x_i - mean)^T = swxx - 2 sw outer(mean, mean)
    + sw outer(mean, mean) = swxx - sw outer(mean, mean).
    """
    n = block.sw
    mean = block.swx / n
    var = np.diag(block.swxx) / n - mean**2
    sd = np.sqrt(np.maximum(var, 0.0))
    sd[sd == 0] = 1.0
    centered_gram = block.swxx - n * np.outer(mean, mean)
    return StandardizedSystem(
        g_std=centered_gram / np.outer(sd, sd),
        b_std=(block.swxy - mean * block.swy) / sd,
        mean=mean,
        sd=sd,
        mean_y=block.swy / n,
    )


def solve(system: StandardizedSystem, alpha: float) -> np.ndarray:
    """Solve (g_std + alpha I) beta_std = b_std for beta_std."""
    p = system.g_std.shape[0]
    return np.linalg.solve(system.g_std + alpha * np.eye(p), system.b_std)


def heldout_r2(
    held: ChromRidgeBlock,
    beta_std: np.ndarray,
    train_mean: np.ndarray,
    train_sd: np.ndarray,
    train_mean_y: float,
) -> float:
    """Weighted R^2 on a held-out block of a model fit on a training set.

    The training model predicts pred_i = c + z_i . beta_std, with c = train_mean_y and
    z_ij = (x_ij - m_j) / s_j using the TRAINING mean m and sd s. Over the held-out
    SNPs, from the held-out block's sufficient statistics alone:

        SS_res = sum_i w_i (y_i - c - z_i . beta_std)^2
               = sum_i w_i (y_i - c)^2 - 2 beta_std . r + beta_std^T Z beta_std
        sum_i w_i (y_i - c)^2 = swyy - 2 c swy + sw c^2
        r_j  = sum_i w_i z_ij (y_i - c)
             = (swxy_j - c swx_j - m_j swy + sw c m_j) / s_j
        Z_jk = sum_i w_i z_ij z_ik
             = (swxx - outer(m, swx) - outer(swx, m) + sw outer(m, m))[j,k] / (s_j s_k)
        SS_tot = sum_i w_i (y_i - ybar)^2 = swyy - sw ybar^2, with ybar = swy / sw

    and R^2 = 1 - SS_res / SS_tot. SS_res is taken about the training intercept c,
    while SS_tot uses the held-out block's own mean ybar: the R^2 null model is
    "predict the held-out mean".
    """
    n = held.sw
    c = train_mean_y
    m, s = train_mean, train_sd
    ss_res_y = held.swyy - 2.0 * c * held.swy + n * c * c
    r = (held.swxy - c * held.swx - m * held.swy + n * c * m) / s
    zz = (
        held.swxx - np.outer(m, held.swx) - np.outer(held.swx, m) + n * np.outer(m, m)
    ) / np.outer(s, s)
    ss_res = ss_res_y - 2.0 * float(beta_std @ r) + float(beta_std @ zz @ beta_std)
    ybar = held.swy / n
    ss_tot = held.swyy - n * ybar * ybar
    return 1.0 - ss_res / ss_tot


def select_alpha_loco(
    blocks_by_chrom: Mapping[int, ChromRidgeBlock], alphas: Sequence[float]
) -> AlphaSelection:
    """Choose alpha maximizing the mean held-out weighted R^2 over
    leave-one-chromosome-out folds of the given blocks."""
    assert len(blocks_by_chrom) >= 2, "LOCO needs at least two chromosome blocks"
    assert len(alphas) > 0, "need at least one candidate alpha"
    assert len(set(alphas)) == len(alphas), "candidate alphas must be distinct"
    chroms = sorted(blocks_by_chrom)
    r2_per_chrom_by_alpha: dict[float, dict[int, float]] = {}
    mean_r2_by_alpha: dict[float, float] = {}
    for alpha in alphas:
        r2s: dict[int, float] = {}
        for held in chroms:
            train = combine([blocks_by_chrom[c] for c in chroms if c != held])
            system = standardized_system(train)
            beta_std = solve(system, alpha)
            r2s[held] = heldout_r2(
                blocks_by_chrom[held],
                beta_std=beta_std,
                train_mean=system.mean,
                train_sd=system.sd,
                train_mean_y=system.mean_y,
            )
        r2_per_chrom_by_alpha[alpha] = r2s
        mean_r2_by_alpha[alpha] = float(np.mean(list(r2s.values())))
    # max keeps the first of tied maxima, so ties go to the earlier candidate.
    best_alpha = max(alphas, key=lambda a: mean_r2_by_alpha[a])
    return AlphaSelection(
        alpha=best_alpha,
        mean_r2=mean_r2_by_alpha[best_alpha],
        r2_per_chrom=r2_per_chrom_by_alpha[best_alpha],
        mean_r2_by_alpha=mean_r2_by_alpha,
    )


def fit(block: ChromRidgeBlock, alpha: float) -> RidgeFit:
    """Fit ridge with penalty alpha on the block and return raw-scale coefficients."""
    system = standardized_system(block)
    beta_std = solve(system, alpha)
    beta_raw = beta_std / system.sd
    intercept = float(system.mean_y - beta_raw @ system.mean)
    return RidgeFit(beta_raw=beta_raw, beta_std=beta_std, intercept=intercept)
