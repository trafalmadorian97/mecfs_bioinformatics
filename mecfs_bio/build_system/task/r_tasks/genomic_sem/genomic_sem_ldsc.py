"""
Pure-Python re-implementation of GenomicSEM's multivariable LD-score
regression (``GenomicSEM::ldsc``).

This produces the same ``covstruc`` that the R function returns:

    S        (k, k)        genetic covariance matrix (liability scale)
    V        (n_V, n_V)    sampling covariance of vech(S), n_V = k(k+1)/2
    I        (k, k)        LDSC intercept / cross-trait intercept matrix
    N        (n_V,)        per-element sample sizes (N_bar for h2, sqrt(N1 N2)
                           for covariances)
    m        float         total number of reference SNPs (sum of M_5_50)
    S_Stand / V_Stand      the standardised (genetic-correlation) versions

The algorithm is a line-by-line port of the installed ``ldsc`` body, so the
conventions match exactly. Notable details that must be preserved for
agreement with R:

- **Weights.** The univariate (heritability) regression weights each SNP by
  ``sqrt(het.w * oc.w)`` normalised to sum 1, where ``het.w`` uses an
  aggregate-h2 estimate and ``oc.w = 1/wLD``. The covariance regression uses
  the *average* of the two traits' initial weights (``weights_cov``). In the
  R source the covariance branch references ``merged$weights`` which does not
  exist and silently partial-matches ``merged$weights_cov`` — so both the
  design and response are weighted by ``weights_cov``.

- **Block jackknife.** ``n_blocks`` contiguous (genome-ordered) blocks. The
  full regression is the sum of per-block X'y and X'X; each leave-one-block-out
  fit gives a delete value, and pseudo-values are
  ``n_blocks * reg - (n_blocks - 1) * delete``. ``V.hold`` stores the
  pseudo-values of the *raw* LD slope (before dividing by N_bar).

- **V scaling.** ``v_out[a,b] = cov(V.hold)[a,b] / (N_a N_b n_blocks / m^2)``,
  which converts the raw-slope jackknife covariance to the covariance of the
  observed-scale h2/cov estimates. Liability scaling is then applied to both
  S and V via the outer product of ``sqrt(Liab.S)``.

- **vech order.** Pairs are enumerated ``(0,0), (0,1), ..., (0,k-1), (1,1),
  ...`` (row-major upper triangle). For a symmetric matrix this indexes the
  same values as R's column-major ``lowerTriangle``, and matches the V_LD
  ordering expected by the GWAS kernels (e.g. (S00, S01, S11) for k=2).
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from attrs import frozen
from scipy.stats import norm

from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_config import (
    LDSC_BP_COL,
    LDSC_CHR_COL,
    LDSC_L2_COL,
    MUNGE_A1_COL,
    MUNGE_N_COL,
    MUNGE_SNP_COL,
    MUNGE_Z_COL,
)

# Derived columns local to this module (not part of the munge/LD-score schemas):
# wLD is the regression-weight LD score, and the cross-trait covariance merge
# suffixes each trait's columns with ".x" (first trait) / ".y" (second).
_WLD_COL = "wLD"
_X_SUFFIX = ".x"
_Y_SUFFIX = ".y"
_Z_X_COL = f"{MUNGE_Z_COL}{_X_SUFFIX}"
_Z_Y_COL = f"{MUNGE_Z_COL}{_Y_SUFFIX}"
_N_X_COL = f"{MUNGE_N_COL}{_X_SUFFIX}"
_N_Y_COL = f"{MUNGE_N_COL}{_Y_SUFFIX}"
_A1_X_COL = f"{MUNGE_A1_COL}{_X_SUFFIX}"
_A1_Y_COL = f"{MUNGE_A1_COL}{_Y_SUFFIX}"


@frozen
class LDSCResult:
    """Output of :func:`run_ldsc`, mirroring GenomicSEM's covstruc list."""

    S: np.ndarray  # (k, k) genetic covariance (liability scale)
    V: np.ndarray  # (n_V, n_V) sampling covariance of vech(S)
    I: np.ndarray  # (k, k) intercepts
    N: np.ndarray  # (n_V,) sample sizes
    m: float  # total reference SNP count
    S_Stand: np.ndarray | None  # (k, k) genetic correlation matrix
    V_Stand: np.ndarray | None  # (n_V, n_V) sampling cov of vech(S_Stand)


@frozen
class _PairEstimate:
    """Result of one heritability or covariance LDSC regression.

    NOTE:
        - reg_tot is scaled slope, so that we estimate h^2 or genetic covariance
        - pseudo_coef are unscaled.  So they are not on the same scale as reg_tot.
    """

    reg_tot: float  # observed-scale h2 (diagonal) or genetic covariance
    intercept: float  # LDSC (or cross-trait) intercept
    n_bar: float  # N_bar (h2) or sqrt(N1 N2) (covariance)
    pseudo_coef: np.ndarray  # (n_blocks,) pseudo-values of the raw LD slope


def _liability_conversion_factor(
    sample_prev: float | None, population_prev: float | None
) -> float:
    """
    Observed -> liability scale conversion factor for a binary trait. Returns
    1.0 (no conversion) when either prevalence is missing.


    See equation 23 in Lee, S. H., Wray, N. R., Goddard, M. E., & Visscher, P. M. (2011). Estimating missing heritability for disease from genome-wide association studies.
    The American Journal of Human Genetics, 88(3), 294-305. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3059431/
    """
    if (
        sample_prev is None
        or population_prev is None
        or np.isnan(sample_prev)
        or np.isnan(population_prev)
    ):
        return 1.0
    pop = population_prev
    samp = sample_prev
    z = float(norm.pdf(norm.ppf(1.0 - pop)))
    return (pop**2 * (1.0 - pop) ** 2) / (samp * (1.0 - samp) * z**2)


def block_bounds(n_snps: int, n_blocks: int) -> list[tuple[int, int]]:
    """
    Contiguous block boundaries matching R's
    ``floor(seq(1, n, length.out = n_blocks + 1))`` scheme. Returns a list of
    ``(start, end)`` half-open Python slice indices.
    """
    assert n_snps >= n_blocks, f"need >= {n_blocks} SNPs; got {n_snps}"
    # 1-indexed block starts, length n_blocks + 1.
    sf = np.floor(np.linspace(1, n_snps, n_blocks + 1)).astype(int)
    bounds: list[tuple[int, int]] = []
    for p in range(n_blocks):
        start_incl_1 = int(sf[p])
        end_incl_1 = int(sf[p + 1]) - 1 if p < n_blocks - 1 else n_snps
        bounds.append((start_incl_1 - 1, end_incl_1))
    return bounds


def _regress_jackknife(
    weighted_ld: np.ndarray, weighted_chi: np.ndarray, n_blocks: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Weighted LD-score regression with a contiguous block jackknife.

    weighted_ld  : (n, 2) design (weighted [L2, intercept])
    weighted_chi : (n,) weighted response
    Returns (reg, pseudo) where reg is (2,) = [slope, intercept] and pseudo is
    (n_blocks, 2) of pseudo-values.


    NOTE:
    reg: Regression coefficient alpha computed via normal equations (X^TX)\alpha = X^TY
    pseudo: Vector of Tukey Jackknife Pseudo-values.  Used to estimate variances and covariances.
    """
    assert weighted_ld.ndim == 2 and weighted_ld.shape[1] == 2
    assert weighted_chi.shape == (weighted_ld.shape[0],)
    bounds = block_bounds(weighted_ld.shape[0], n_blocks)

    block_xty = np.empty((n_blocks, 2))
    block_xtx = np.empty((n_blocks, 2, 2))
    for i, (a, b) in enumerate(bounds):
        wx = weighted_ld[a:b]
        wy = weighted_chi[a:b]
        block_xty[i] = wx.T @ wy
        block_xtx[i] = wx.T @ wx

    xty = block_xty.sum(
        axis=0
    )  # This amounts to computing X^TX from the normal equations through block matrix multiplication
    xtx = block_xtx.sum(
        axis=0
    )  # This amounts to computing X^TY from the normal equations through block matrix multiplication
    reg = np.linalg.solve(xtx, xty)

    pseudo = np.empty((n_blocks, 2))
    for i in range(n_blocks):
        delete = np.linalg.solve(xtx - block_xtx[i], xty - block_xty[i])
        pseudo[i] = n_blocks * reg - (n_blocks - 1) * delete
    return reg, pseudo


def _aggregate_h2(chi: np.ndarray, ld: np.ndarray, n: np.ndarray, m: float) -> float:
    """Aggregate h2 estimate used to build heteroskedasticity weights, clamped
    to [0, 1] (R's tot.agg).


    NOTE:
    Recall that the core LDSC equation is

    E chi^2 = h^2*l*N/m +1

    The h^2 estimate computed by this function comes from taking the mean of the key quantities on both sides and solving
    for h^2.


    """
    tot_agg = (m * (float(np.mean(chi)) - 1.0)) / float(np.mean(ld * n))
    return float(min(max(tot_agg, 0.0), 1.0))


def _het_oc_initial_weight(
    tot_agg: float, ld_raw: np.ndarray, wld_raw: np.ndarray, n: np.ndarray, m: float
) -> np.ndarray:
    """sqrt(het.w * oc.w) per SNP (R's initial.w).



    NOTE:
    The purpose of this function is to compute regression weights for use in the LDSC regression.

    There are two reasons for computing these regression weights:
    1.  Heteroskedasticity: The signal from different SNPs will have different variances.  We should downweight SNPs with
        high variances.

    2.  Overcounting: because neighbouring SNPs are in linkage disequilibrium, they to some extent carry redundant information.
        We use weighting to avoid overcounting this redundant information.


    The het_w weights computed by this function account for heteroskedasticity.  The oc_w weights account for overcounting.

    The computation of the het_w weights takes advantage of the fact that in the LDSC model, the chi^2 statistics are distributed
    according to a standard chi^2 distributed scaled by factor = h^2*m/n*l_j +1.  A standard chi^2 random variable has
    variance of twice its mean.  Scaling the random variable by a factor scales its variance by factor**2.

    """
    ld = np.maximum(ld_raw, 1.0)
    w_ld = np.maximum(wld_raw, 1.0)
    c = tot_agg * n / m
    het_w = 1.0 / (2.0 * (1.0 + c * ld) ** 2)
    oc_w = 1.0 / w_ld
    return np.sqrt(het_w * oc_w)


def estimate_h2(
    *,
    chi: np.ndarray,
    ld_raw: np.ndarray,
    wld_raw: np.ndarray,
    n: np.ndarray,
    m: float,
    n_blocks: int,
) -> _PairEstimate:
    tot_agg = _aggregate_h2(chi, ld_raw, n, m)
    initial_w = _het_oc_initial_weight(tot_agg, ld_raw, wld_raw, n, m)
    weights = initial_w / initial_w.sum()
    n_bar = float(np.mean(n))

    weighted_ld = np.column_stack([ld_raw, np.ones_like(ld_raw)]) * weights[:, None]
    weighted_chi = chi * weights
    reg, pseudo = _regress_jackknife(weighted_ld, weighted_chi, n_blocks)

    intercept = float(reg[1])
    reg_tot = float(reg[0] / n_bar * m)
    return _PairEstimate(
        reg_tot=reg_tot, intercept=intercept, n_bar=n_bar, pseudo_coef=pseudo[:, 0]
    )


def estimate_cov(
    *,
    zz: np.ndarray,
    chi1: np.ndarray,
    chi2: np.ndarray,
    ld_raw: np.ndarray,
    wld_raw: np.ndarray,
    n_x: np.ndarray,
    n_y: np.ndarray,
    m: float,
    n_blocks: int,
) -> _PairEstimate:
    tot_agg = _aggregate_h2(chi1, ld_raw, n_x, m)
    tot_agg2 = _aggregate_h2(chi2, ld_raw, n_y, m)
    initial_w = _het_oc_initial_weight(tot_agg, ld_raw, wld_raw, n_x, m)
    initial_w2 = _het_oc_initial_weight(tot_agg2, ld_raw, wld_raw, n_y, m)
    # R references merged$weights here, which partial-matches weights_cov.
    weights_cov = (initial_w + initial_w2) / (initial_w + initial_w2).sum()
    n_bar = float(np.sqrt(np.mean(n_x) * np.mean(n_y)))

    weighted_ld = np.column_stack([ld_raw, np.ones_like(ld_raw)]) * weights_cov[:, None]
    weighted_chi = zz * weights_cov
    reg, pseudo = _regress_jackknife(weighted_ld, weighted_chi, n_blocks)

    intercept = float(reg[1])
    reg_tot = float(reg[0] / n_bar * m)
    return _PairEstimate(
        reg_tot=reg_tot, intercept=intercept, n_bar=n_bar, pseudo_coef=pseudo[:, 0]
    )


def _read_ld_scores(ld_dir: Path, n_chrom: int) -> tuple[pd.DataFrame, float]:
    """
    Read ``<chr>.l2.ldscore.gz`` (SNP, CHR, BP, L2) and ``<chr>.l2.M_5_50``
    for chromosomes 1..n_chrom. Returns (ld_df, m_total).

    Point to an LD directory in standard format
    - Reads all dataframes and concatenates them.  Return the concatenated result
    - Reads all the M_5_50 files, which store the number of SNPs with minor allele frequencies between 5 and 50 used to compute the LD scores
      Sums these values across all chromosomes, to return an M total.
    """
    ld_frames = []
    m_total = 0.0
    for chrom in range(1, n_chrom + 1):
        score = pd.read_csv(ld_dir / f"{chrom}.l2.ldscore.gz", sep="\t")
        ld_frames.append(score[[LDSC_CHR_COL, MUNGE_SNP_COL, LDSC_BP_COL, LDSC_L2_COL]])
        m_val = pd.read_csv(ld_dir / f"{chrom}.l2.M_5_50", header=None)
        m_total += float(np.asarray(m_val).sum())
    ld_df = pd.concat(ld_frames, ignore_index=True)
    return ld_df, m_total


def _read_munged(path: Path) -> pd.DataFrame:
    """Read a munged ``.sumstats(.gz)`` file, keeping SNP, N, Z, A1."""
    df = pd.read_csv(path, sep="\t")
    df = df[[MUNGE_SNP_COL, MUNGE_N_COL, MUNGE_Z_COL, MUNGE_A1_COL]].dropna()
    return df


def _merge_trait_with_ld_and_sort(
    trait: pd.DataFrame, ld_df: pd.DataFrame, chisq_max: float | None
) -> pd.DataFrame:
    """
    Merge one trait's sumstats with LD scores, order by (CHR, BP), and apply
    the chi^2 filter (per-trait threshold when chisq_max is None).
    """
    merged = trait.merge(ld_df, on=MUNGE_SNP_COL, how="inner", sort=False)
    merged = merged.sort_values([LDSC_CHR_COL, LDSC_BP_COL], kind="stable").reset_index(
        drop=True
    )
    threshold = (
        chisq_max
        if chisq_max is not None
        else max(0.001 * float(merged[MUNGE_N_COL].max()), 80.0)
    )
    keep = merged[MUNGE_Z_COL].to_numpy() ** 2 <= threshold
    return merged.loc[keep].reset_index(drop=True)


def run_ldsc(
    *,
    munged_paths: Sequence[Path],
    ld_dir: Path,
    sample_prev: Sequence[float | None],
    population_prev: Sequence[float | None],
    n_blocks: int = 200,
    n_chrom: int = 22,
    chisq_max: float | None = None,
    stand: bool = True,
) -> LDSCResult:
    """
    Multivariable LD-score regression, faithful to ``GenomicSEM::ldsc``.

    munged_paths    : per-trait munged sumstats files (order defines trait order)
    ld_dir          : directory of ``<chr>.l2.ldscore.gz`` / ``<chr>.l2.M_5_50``
    sample_prev     : per-trait sample prevalence (None for continuous traits)
    population_prev : per-trait population prevalence (None for continuous)
    """
    k = len(munged_paths)
    assert k >= 2, "ldsc requires >= 2 traits"
    assert len(sample_prev) == len(population_prev) == k

    ld_df, m_total = _read_ld_scores(ld_dir, n_chrom)
    # sep_weights = FALSE: weights LD score equals the regression LD score.
    ld_df = ld_df.assign(**{_WLD_COL: ld_df[LDSC_L2_COL]})

    trait_frames = [
        _merge_trait_with_ld_and_sort(_read_munged(Path(p)), ld_df, chisq_max)
        for p in munged_paths
    ]

    liab = np.array(
        [
            _liability_conversion_factor(sample_prev[j], population_prev[j])
            for j in range(k)
        ]
    )

    # Enumerate vech pairs in row-major upper-triangle order.
    pairs = [(j, kk) for j in range(k) for kk in range(j, k)]
    n_v = len(pairs)

    cov = np.full((k, k), np.nan)
    intercepts = np.full((k, k), np.nan)
    n_vec = np.empty(n_v)
    v_hold = np.empty((n_blocks, n_v))

    for s, (j, kk) in enumerate(pairs):
        if j == kk:
            yj = trait_frames[j]
            est = estimate_h2(
                chi=yj[MUNGE_Z_COL].to_numpy() ** 2,
                ld_raw=yj[LDSC_L2_COL].to_numpy(),
                wld_raw=yj[_WLD_COL].to_numpy(),
                n=yj[MUNGE_N_COL].to_numpy(dtype=float),
                m=m_total,
                n_blocks=n_blocks,
            )
            cov[j, j] = est.reg_tot
            intercepts[j, j] = est.intercept
        else:
            merged = trait_frames[j].merge(
                trait_frames[kk][
                    [MUNGE_SNP_COL, MUNGE_N_COL, MUNGE_Z_COL, MUNGE_A1_COL]
                ],
                on=MUNGE_SNP_COL,
                how="inner",
                sort=False,
                suffixes=(_X_SUFFIX, _Y_SUFFIX),
            )
            merged = merged.sort_values(
                [LDSC_CHR_COL, LDSC_BP_COL], kind="stable"
            ).reset_index(drop=True)
            z_x = merged[_Z_X_COL].to_numpy()
            z_y = merged[_Z_Y_COL].to_numpy()
            # Flip trait-j Z where the reference alleles disagree.
            same_allele = (merged[_A1_Y_COL] == merged[_A1_X_COL]).to_numpy()
            z_x = np.where(same_allele, z_x, -z_x)
            est = estimate_cov(
                zz=z_y * z_x,
                chi1=z_x**2,
                chi2=z_y**2,
                ld_raw=merged[LDSC_L2_COL].to_numpy(),
                wld_raw=merged[_WLD_COL].to_numpy(),
                n_x=merged[_N_X_COL].to_numpy(dtype=float),
                n_y=merged[_N_Y_COL].to_numpy(dtype=float),
                m=m_total,
                n_blocks=n_blocks,
            )
            cov[j, kk] = cov[kk, j] = est.reg_tot
            intercepts[j, kk] = intercepts[kk, j] = est.intercept

        n_vec[s] = est.n_bar
        v_hold[:, s] = est.pseudo_coef

    # Raw-slope jackknife covariance -> observed-scale h2/cov sampling cov.
    denom = np.outer(
        n_vec * (np.sqrt(n_blocks) / m_total), n_vec * (np.sqrt(n_blocks) / m_total)
    )  # Note: this denominator combines two re-scaling. One prescribed by Jackknife, and one to convert
    # raw regression slope to heritability/ covariance units
    v_out = np.cov(v_hold, rowvar=False) / denom

    # Liability scaling of S and V.
    sqrt_liab = np.sqrt(liab)
    ratio = np.outer(sqrt_liab, sqrt_liab)
    S = cov * ratio
    scale_o = np.array([ratio[j, kk] for (j, kk) in pairs])
    V = v_out * np.outer(scale_o, scale_o)

    S_stand: np.ndarray | None = None
    V_stand: np.ndarray | None = None
    if stand and np.all(np.diag(S) > 0):
        inv_sd = 1.0 / np.sqrt(np.diag(S))
        ratio_stand = np.outer(inv_sd, inv_sd)
        S_stand = S * ratio_stand
        scale_o_stand = np.array([ratio_stand[j, kk] for (j, kk) in pairs])
        V_stand = V * np.outer(scale_o_stand, scale_o_stand)

    return LDSCResult(
        S=S,
        V=V,
        I=intercepts,
        N=n_vec,
        m=m_total,
        S_Stand=S_stand,
        V_Stand=V_stand,
    )
