"""
Batched cross-trait LD-score-regression genetic correlation for the UKB-PPP database.

One external trait is correlated against every protein. Because every protein's z-score is
defined on the SAME shared SNP set (see ppp_ldsc_context), we process K proteins at once and
run the three regressions each genetic correlation needs:

    rg = gcov(trait, protein) / sqrt(h2_trait * h2_protein)

Following the GenomicSEM convention (see genomic_sem_ldsc.run_ldsc), each DIAGONAL heritability
is estimated on that trait's own kept SNP set, and only the OFF-DIAGONAL covariance on the
pair's intersection. So h2_trait -- its point estimate, its per-block leave-one-out delete
values, and n_bar -- depends only on the trait and is computed exactly ONCE
(estimate_trait_context) and reused for every protein. Only h2_protein and gcov are per-protein.

Per-protein exclusions (missing variants, the chi-square filter, and cis variants) are encoded
as ZERO regression weights over the shared row set rather than by dropping rows, so all proteins
keep the same contiguous jackknife blocks -- exactly as in batched_ldsc_h2.

Weighting (the easy bug): the univariate h2 weight is omega = het * oc (the PRODUCT), because
estimate_h2 multiplies BOTH design and response by sqrt(het * oc). The covariance weight is
omega_cov = oc * (sqrt(het_trait) + sqrt(het_protein))^2, from estimate_cov's
weights_cov = (initial_w_trait + initial_w_protein) with initial_w = sqrt(het * oc); the
sum-to-1 normalization is a common scalar that cancels in the 2x2 slope solve. The heteroskedasticity
aggregate-h2 that shapes het is computed over the SAME set the regression runs on: own kept set
for a diagonal h2, the pair intersection for the covariance.

The rg standard error is a delete-block jackknife on the ratio: leave-one-block-out slopes for
gcov and h2_protein are combined with the once-computed h2_trait deletes at the same block index
to form per-block rg values, and rg_se is their jackknife spread. This is the classic ldsc rg SE.
"""

from __future__ import annotations

import numpy as np
from attrs import frozen
from scipy.stats import norm

from mecfs_bio.build_system.task.ppp_ldsc.batched_ldsc_h2 import DEFAULT_N_BLOCKS
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_ldsc import (
    block_bounds,
    estimate_cov,
    estimate_h2,
)

_CHISQ_FILTER_FLOOR = 80.0
_CHISQ_FILTER_N_COEF = 0.001


def _chisq_threshold(n: np.ndarray | float) -> np.ndarray:
    """Per-trait chi-square filter max(0.001 * N, 80), matching GenomicSEM. Returns an array
    for array input (subscriptable per protein) and a 0-d array for scalar input."""
    return np.asarray(
        np.maximum(
            _CHISQ_FILTER_N_COEF * np.asarray(n, dtype=float), _CHISQ_FILTER_FLOOR
        )
    )


def _solve(sums: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    r"""Solve the 2x2 weighted normal equations for design [ld, 1] from the packed block sums
       (a, b, d, e, f) = (sum w*ld^2, sum w*ld, sum w, sum w*ld*y, sum w*y). Returns (slope,
       intercept), broadcasting over any leading axes.

       NOTES:

       This solves the 2x2 weighted normal equations.
       We have:

       The weighted normal equations are

       A^TWy = A^TWAb

       In our case:

       A = [L  1] \in \R^{nx2}
       A^TWY = [L^T]   W   y
               [1^T]
             = [L^TWy]
               [1^TWy]

             =:[e]     \in R^2
               [f]


       A^TWA  =  [L^T] W   [ L  1]
                 [1^T]
              =  [ L^TWL     L^TW1]
                [ 1^TWL     1^Twq]
              =: [a   b]    \in R^{2x2}
                 [c   d]



        (A^TWA)^{-1} = 1/(ad-b^2)    [d    -b]
                                     [-b    a]


        (A^TWA)^{-1}A^TWY = 1/(ad-b^2)  [de-bf]
                                        [af-be]

    Which is the solution of the normal equations
    """
    a, b, d, e, f = (sums[..., j] for j in range(5))
    det = a * d - b * b
    slope = (d * e - b * f) / det
    intercept = (a * f - b * e) / det
    return slope, intercept


def _blocked_wls(
    w: np.ndarray, y: np.ndarray, ld: np.ndarray, n_blocks: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Weighted LD-score regression of y on [ld, 1] with a contiguous block jackknife, for a
    stack of K columns sharing the row set. w, y are (S, K); ld is (S,). Returns
    (slope_full (K,), intercept_full (K,), slope_delete (n_blocks, K)) where slope_delete[i] is
    the leave-block-i-out slope."""
    s, k = w.shape
    ld_col = ld[:, None]
    ld2 = (ld * ld)[:, None]
    blk = np.empty((n_blocks, k, 5))
    for i, (lo, hi) in enumerate(block_bounds(s, n_blocks)):
        wb = w[lo:hi]
        yb = y[lo:hi]
        blk[i, :, 0] = (wb * ld2[lo:hi]).sum(0)
        blk[i, :, 1] = (wb * ld_col[lo:hi]).sum(0)
        blk[i, :, 2] = wb.sum(0)
        blk[i, :, 3] = (wb * ld_col[lo:hi] * yb).sum(0)
        blk[i, :, 4] = (wb * yb).sum(0)
    tot = blk.sum(0)
    slope_full, intercept_full = _solve(tot)
    slope_delete, _ = _solve(tot[None, :, :] - blk)
    return slope_full, intercept_full, slope_delete


def _aggregate_h2(
    chi: np.ndarray, ld: np.ndarray, n: np.ndarray, keep: np.ndarray, m: float
) -> np.ndarray:
    """Aggregate-h2 per column over the kept SNPs (GenomicSEM tot.agg), clamped to [0, 1]. chi,
    n, keep are (S, K); ld is (S,). Uses raw ld (not clamped) for mean(ld*N), matching
    genomic_sem_ldsc._aggregate_h2.

    Note: idea is to compute a rough estimate of heritability that can be used for variant weighting
    """
    cnt = np.maximum(keep.sum(0).astype(float), 1.0)
    mean_chi = np.where(keep, chi, 0.0).sum(0) / cnt
    mean_ldn = np.where(keep, ld[:, None] * n, 0.0).sum(0) / cnt
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.clip(m * (mean_chi - 1.0) / mean_ldn, 0.0, 1.0)


def _sqrt_het(
    tot_agg: np.ndarray, ldm: np.ndarray, n: np.ndarray, m: float
) -> np.ndarray:
    """Per-SNP sqrt of the heteroskedasticity weight het = 1/(2*(1 + c*ldm)^2) with
    c = tot_agg*N/m, mirroring genomic_sem_ldsc._het_oc_initial_weight. tot_agg is (K,), ldm is
    (S,), n is (S, K); returns (S, K)."""
    c = tot_agg[None, :] * n / m
    het = 1.0 / (2.0 * (1.0 + c * ldm[:, None]) ** 2)
    return np.sqrt(het)


def _column_mean_over_keep(
    value: np.ndarray, keep: np.ndarray, cnt: np.ndarray
) -> np.ndarray:
    """Mean of value (S, K) over the kept SNPs per column; NaN where a column keeps nothing."""
    return np.where(
        cnt > 0, np.where(keep, value, 0.0).sum(0) / np.maximum(cnt, 1.0), np.nan
    )


@frozen
class TraitLdscContext:
    """The trait-only LDSC quantities, computed once and reused for every protein.

    z / keep / n are parallel (S,) arrays over the shared context SNP set: z is the signed
    trait z-score (index-effect-allele oriented, zeroed off the kept set), keep marks the trait's
    kept SNPs (present and passing the chi-square filter), and n is the per-SNP trait sample size.
    h2 / n_bar are the trait heritability and its N_bar; h2_delete holds the (n_blocks,) per-block
    leave-one-out h2_trait values used by the rg jackknife.
    """

    z: np.ndarray
    keep: np.ndarray
    n: np.ndarray
    h2: float
    n_bar: float
    h2_delete: np.ndarray

    def __attrs_post_init__(self) -> None:
        s = self.z.shape[0]
        for name, arr, kind in (
            ("z", self.z, "f"),
            ("keep", self.keep, "b"),
            ("n", self.n, "f"),
        ):
            assert arr.ndim == 1, f"{name} must be 1-D, got shape {arr.shape}"
            assert arr.shape[0] == s, f"{name} has length {arr.shape[0]}, expected {s}"
            assert arr.dtype.kind == kind, (
                f"{name} must be {kind!r}-kind, got dtype {arr.dtype}"
            )
        assert self.h2_delete.ndim == 1, "h2_delete must be 1-D"


@frozen
class BatchedRgResult:
    """Per-protein cross-trait rg outputs; all arrays are length K (protein order preserved).
    h2_trait and n_bar_trait are trait-level scalars shared by every protein."""

    rg: np.ndarray
    rg_se: np.ndarray
    rg_z: np.ndarray
    rg_p: np.ndarray
    gcov: np.ndarray
    gcov_intercept: np.ndarray
    h2_protein: np.ndarray
    n_snps: np.ndarray  # SNPs entering the covariance (trait-protein intersection)
    n_bar_protein: np.ndarray
    h2_trait: float
    n_bar_trait: float

    def __attrs_post_init__(self) -> None:
        k = self.rg.shape[0]
        for name, arr, kind in (
            ("rg", self.rg, "f"),
            ("rg_se", self.rg_se, "f"),
            ("rg_z", self.rg_z, "f"),
            ("rg_p", self.rg_p, "f"),
            ("gcov", self.gcov, "f"),
            ("gcov_intercept", self.gcov_intercept, "f"),
            ("h2_protein", self.h2_protein, "f"),
            ("n_snps", self.n_snps, "i"),
            ("n_bar_protein", self.n_bar_protein, "f"),
        ):
            assert arr.ndim == 1, f"{name} must be 1-D, got shape {arr.shape}"
            assert arr.shape[0] == k, f"{name} has length {arr.shape[0]}, expected {k}"
            assert arr.dtype.kind == kind, (
                f"{name} must be {kind!r}-kind, got dtype {arr.dtype}"
            )


def estimate_trait_context(
    z_trait: np.ndarray,
    n_trait: np.ndarray,
    ld: np.ndarray,
    m: float,
    *,
    n_blocks: int = DEFAULT_N_BLOCKS,
) -> TraitLdscContext:
    """Estimate the trait heritability once over the shared SNP set and package the trait-only
    quantities the batched rg kernel reuses for every protein.

    z_trait: (S,) signed trait z-score aligned to the context (index) effect allele, NaN where
        the trait lacks the variant. n_trait: (S,) trait sample size (NaN where absent). ld: (S,)
        LD score, genome-sorted for contiguous blocks. m: total reference-SNP count.
    """
    finite = np.isfinite(z_trait)
    n_col = np.where(finite, n_trait, 0.0)
    chi = np.where(finite, z_trait * z_trait, 0.0)
    # The per-trait chi-square threshold uses the trait's max N (as in genomic_sem_ldsc).
    threshold = _chisq_threshold(float(np.nanmax(n_trait)))
    keep = finite & (chi <= threshold)
    z = np.where(keep, z_trait, 0.0)

    keep_col = keep[:, None]
    n_sk = n_col[:, None]
    ldm = np.maximum(ld, 1.0)
    oc = 1.0 / ldm
    tot_agg = _aggregate_h2(chi[:, None], ld, n_sk, keep_col, m)
    sqrt_het = _sqrt_het(tot_agg, ldm, n_sk, m)
    w = sqrt_het**2 * oc[:, None] * keep_col
    slope, _intercept, slope_delete = _blocked_wls(w, chi[:, None], ld, n_blocks)

    cnt = float(keep.sum())
    n_bar = float(np.where(keep, n_col, 0.0).sum() / max(cnt, 1.0))
    h2 = float(slope[0] / n_bar * m)
    h2_delete = slope_delete[:, 0] / n_bar * m
    return TraitLdscContext(
        z=z.astype(float),
        keep=keep,
        n=n_col.astype(float),
        h2=h2,
        n_bar=n_bar,
        h2_delete=h2_delete,
    )


def batched_rg(
    trait_ctx: TraitLdscContext,
    z_protein: np.ndarray,
    ld: np.ndarray,
    n_protein: np.ndarray,
    m: float,
    *,
    n_blocks: int = DEFAULT_N_BLOCKS,
    exclude: np.ndarray | None = None,
) -> BatchedRgResult:
    """Estimate the trait-vs-protein genetic correlation for K proteins sharing the SNP set.

    trait_ctx: the once-computed trait LDSC context (estimate_trait_context).
    z_protein: (S, K) signed protein z-score (= BETA/SE, index-effect-allele oriented), NaN where
        a protein lacks the variant.
    ld: (S,) LD score, genome-sorted. n_protein: (K,) per-protein sample size. m: total
        reference-SNP count. exclude: optional (S, K) boolean; True drops that variant for that
        protein (e.g. cis).
    """
    s, k = z_protein.shape
    ldm = np.maximum(ld, 1.0)
    oc = 1.0 / ldm

    finite_p = np.isfinite(z_protein)
    z_p = np.where(finite_p, z_protein, 0.0)
    chi_p = z_p * z_p
    n_p = np.broadcast_to(n_protein.astype(float), (s, k))
    keep_p = finite_p & (chi_p <= _chisq_threshold(n_protein)[None, :])
    if exclude is not None:
        keep_p = keep_p & ~exclude

    # Protein heritability: own-set aggregate-h2, weight het*oc (the univariate PRODUCT).
    tot_agg_p = _aggregate_h2(chi_p, ld, n_p, keep_p, m)
    w_p = _sqrt_het(tot_agg_p, ldm, n_p, m) ** 2 * oc[:, None] * keep_p
    slope_p, _intc_p, slope_p_del = _blocked_wls(w_p, chi_p, ld, n_blocks)
    cnt_p = keep_p.sum(0).astype(float)
    n_bar_p = _column_mean_over_keep(n_p, keep_p, cnt_p)
    h2_p = slope_p / n_bar_p * m
    h2_p_del = slope_p_del / n_bar_p[None, :] * m

    # Trait-protein covariance on the intersection: het weights use intersection aggregate-h2.
    keep_cov = trait_ctx.keep[:, None] & keep_p
    chi_t = (trait_ctx.z * trait_ctx.z)[:, None]
    n_t = trait_ctx.n[:, None]
    zz = np.where(keep_cov, trait_ctx.z[:, None] * z_p, 0.0)
    tot_agg_t_cov = _aggregate_h2(
        np.broadcast_to(chi_t, (s, k)), ld, np.broadcast_to(n_t, (s, k)), keep_cov, m
    )
    tot_agg_p_cov = _aggregate_h2(chi_p, ld, n_p, keep_cov, m)
    sqrt_het_t = _sqrt_het(tot_agg_t_cov, ldm, np.broadcast_to(n_t, (s, k)), m)
    sqrt_het_p = _sqrt_het(tot_agg_p_cov, ldm, n_p, m)
    w_cov = oc[:, None] * (sqrt_het_t + sqrt_het_p) ** 2 * keep_cov
    slope_cov, intercept_cov, slope_cov_del = _blocked_wls(w_cov, zz, ld, n_blocks)

    cnt_cov = keep_cov.sum(0).astype(float)
    mean_n_t = _column_mean_over_keep(np.broadcast_to(n_t, (s, k)), keep_cov, cnt_cov)
    mean_n_p = _column_mean_over_keep(n_p, keep_cov, cnt_cov)
    n_bar_cov = np.sqrt(mean_n_t * mean_n_p)
    gcov = slope_cov / n_bar_cov * m
    gcov_del = slope_cov_del / n_bar_cov[None, :] * m

    # rg and its delete-block jackknife SE on the ratio.
    invalid = (h2_p <= 0.0) | (trait_ctx.h2 <= 0.0)
    with np.errstate(invalid="ignore"):
        rg = np.where(invalid, np.nan, gcov / np.sqrt(trait_ctx.h2 * h2_p))
        denom_del = np.sqrt(trait_ctx.h2_delete[:, None] * h2_p_del)
        rg_del = gcov_del / denom_del
        rg_del_mean = rg_del.mean(0)
        rg_se = np.sqrt(
            (n_blocks - 1) / n_blocks * ((rg_del - rg_del_mean[None, :]) ** 2).sum(0)
        )
        rg_se = np.where(invalid, np.nan, rg_se)
        rg_z = rg / rg_se
    rg_p = 2.0 * norm.sf(np.abs(rg_z))

    return BatchedRgResult(
        rg=rg,
        rg_se=rg_se,
        rg_z=rg_z,
        rg_p=rg_p,
        gcov=gcov,
        gcov_intercept=intercept_cov,
        h2_protein=h2_p,
        n_snps=cnt_cov.astype(np.int64),
        n_bar_protein=n_bar_p,
        h2_trait=trait_ctx.h2,
        n_bar_trait=trait_ctx.n_bar,
    )


@frozen
class ExactRgResult:
    rg: float
    gcov: float
    gcov_intercept: float
    h2_trait: float
    h2_protein: float


def exact_rg_single(
    z_trait: np.ndarray,
    n_trait: float,
    z_protein: np.ndarray,
    ld: np.ndarray,
    n_protein: float,
    m: float,
    *,
    n_blocks: int = DEFAULT_N_BLOCKS,
    exclude: np.ndarray | None = None,
) -> ExactRgResult:
    """Reference single-protein genetic correlation via the repo's exact GenomicSEM port: h2 on
    each trait's own kept set and gcov on the intersection (dropping, not zero-weighting, the
    filtered/cis SNPs). This is what the batched kernel's point estimates are validated against.
    z_trait/z_protein/ld are genome-sorted (S,) arrays; n_trait/n_protein are scalar."""
    keep_t = np.isfinite(z_trait) & (z_trait**2 <= _chisq_threshold(n_trait))
    keep_p = np.isfinite(z_protein) & (z_protein**2 <= _chisq_threshold(n_protein))
    if exclude is not None:
        keep_p = keep_p & ~exclude
    keep_cov = keep_t & keep_p

    est_t = estimate_h2(
        chi=z_trait[keep_t] ** 2,
        ld_raw=ld[keep_t],
        wld_raw=ld[keep_t],
        n=np.full(int(keep_t.sum()), n_trait),
        m=m,
        n_blocks=n_blocks,
    )
    est_p = estimate_h2(
        chi=z_protein[keep_p] ** 2,
        ld_raw=ld[keep_p],
        wld_raw=ld[keep_p],
        n=np.full(int(keep_p.sum()), n_protein),
        m=m,
        n_blocks=n_blocks,
    )
    est_c = estimate_cov(
        zz=z_trait[keep_cov] * z_protein[keep_cov],
        chi1=z_trait[keep_cov] ** 2,
        chi2=z_protein[keep_cov] ** 2,
        ld_raw=ld[keep_cov],
        wld_raw=ld[keep_cov],
        n_x=np.full(int(keep_cov.sum()), n_trait),
        n_y=np.full(int(keep_cov.sum()), n_protein),
        m=m,
        n_blocks=n_blocks,
    )
    rg = est_c.reg_tot / np.sqrt(est_t.reg_tot * est_p.reg_tot)
    return ExactRgResult(
        rg=float(rg),
        gcov=est_c.reg_tot,
        gcov_intercept=est_c.intercept,
        h2_trait=est_t.reg_tot,
        h2_protein=est_p.reg_tot,
    )
