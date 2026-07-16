"""
Batched univariate LD-score-regression heritability for the UKB-PPP database.

Because every protein's chi-square is defined on the SAME shared SNP set (see
ppp_ldsc_context), we process K proteins at once: stack chi-square into an (S, K) matrix
and run the weighted LD-score regression plus block jackknife for all K in vectorized
numpy. Per-protein exclusions (missing variants, the chi-square filter, and cis variants)
are encoded as ZERO regression weights rather than by dropping rows, so all proteins keep
the same row set and the same shared contiguous jackknife blocks.

This is validated (experiments/claude/ppp_ldsc/batched_vs_exact_h2_probe.py) to agree with
the repo's exact per-protein GenomicSEM port (_genomic_sem_ldsc._estimate_h2): the h2 point
estimate is machine-identical and the jackknife SE agrees to <=~1.5% (the only difference is
that the exact method cuts the KEPT SNPs into blocks whereas the shared-block method cuts the
FULL set, shifting block membership by a few SNPs).

Weighting note (the one easy bug): the effective WLS weight is omega = het * oc (the
PRODUCT), because _estimate_h2 multiplies BOTH the design and the response by
initial_w = sqrt(het * oc). The sum-to-1 normalization is a common scalar that cancels in
the 2x2 slope solve, so it is omitted here.
"""

from __future__ import annotations

import numpy as np
from attrs import frozen

from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_ldsc import (
    _block_bounds,
    _estimate_h2,
)

DEFAULT_N_BLOCKS = 200
# Median of a 1-df chi-square under the null, for genomic-control lambda_gc.
_CHI2_1DF_MEDIAN = 0.4549364231195724


@frozen
class BatchedH2Result:
    """Per-protein heritability outputs; all arrays are length K (protein order preserved)."""

    h2: np.ndarray
    h2_se: np.ndarray
    intercept: np.ndarray
    mean_chi2: np.ndarray
    lambda_gc: np.ndarray
    n_snps: np.ndarray  # kept SNP count per protein


def _chisq_threshold(n: np.ndarray) -> np.ndarray:
    """Per-protein chi-square filter max(0.001 * N, 80), matching GenomicSEM's per-trait
    threshold (with N constant across a protein's SNPs)."""
    return np.maximum(0.001 * n, 80.0)


def batched_h2(
    chi2: np.ndarray,
    ld: np.ndarray,
    n: np.ndarray,
    m: float,
    *,
    n_blocks: int = DEFAULT_N_BLOCKS,
    exclude: np.ndarray | None = None,
) -> BatchedH2Result:
    """Estimate LDSC heritability for K proteins sharing the SNP set.

    chi2: (S, K) chi-square (= (BETA/SE)^2), NaN where a protein lacks the variant.
    ld: (S,) LD score (== weight LD score), genome-sorted for contiguous blocks.
    n: (K,) per-protein sample size (constant across a protein's SNPs).
    m: total reference-SNP count.
    exclude: optional (S, K) boolean; True drops that variant for that protein (e.g. cis).
    """
    s, k = chi2.shape
    keep = np.isfinite(chi2) & (chi2 <= _chisq_threshold(n)[None, :])
    if exclude is not None:
        keep &= ~exclude
    chi = np.where(keep, chi2, 0.0)

    # Aggregate-h2 per protein over KEPT SNPs (GenomicSEM tot.agg), for the weights.
    cnt = keep.sum(0).astype(float)  # (K,)
    ld_col = ld[:, None]
    mean_chi = chi.sum(0) / cnt
    mean_ldn = (
        np.where(keep, ld_col, 0.0).sum(0) * n
    ) / cnt  # mean(ld)*N == mean(ld*N)
    tot_agg = np.clip(m * (mean_chi - 1.0) / mean_ldn, 0.0, 1.0)  # (K,)

    # Effective WLS weight omega = het * oc (see module docstring), zeroed off the kept set.
    ldm = np.maximum(ld, 1.0)[:, None]  # wLD == LD here
    c = (tot_agg * n / m)[None, :]
    het = 1.0 / (2.0 * (1.0 + c * ldm) ** 2)
    oc = 1.0 / ldm
    w = het * oc * keep  # (S, K)

    # Per-block weighted sums for the 2x2 normal equations of design [ld, 1]:
    #   a = sum w*ld^2, b = sum w*ld, d = sum w ; e = sum w*ld*chi, f = sum w*chi.
    ld2 = (ld * ld)[:, None]
    blk = np.empty((n_blocks, k, 5))
    for i, (lo, hi) in enumerate(_block_bounds(s, n_blocks)):
        wb = w[lo:hi]
        cb = chi[lo:hi]
        blk[i, :, 0] = (wb * ld2[lo:hi]).sum(0)
        blk[i, :, 1] = (wb * ld_col[lo:hi]).sum(0)
        blk[i, :, 2] = wb.sum(0)
        blk[i, :, 3] = (wb * ld_col[lo:hi] * cb).sum(0)
        blk[i, :, 4] = (wb * cb).sum(0)
    tot = blk.sum(0)  # (K, 5)

    def _solve(sums: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        a, b, d, e, f = (sums[..., j] for j in range(5))
        det = a * d - b * b
        slope = (d * e - b * f) / det
        intercept = (a * f - b * e) / det
        return slope, intercept

    slope_full, intercept_full = _solve(tot)

    # Leave-one-block-out slopes -> pseudo-values of the raw slope -> jackknife SE.
    slope_loo, _ = _solve(tot[None, :, :] - blk)  # (n_blocks, K)
    pseudo = n_blocks * slope_full[None, :] - (n_blocks - 1) * slope_loo
    denom = (n * np.sqrt(n_blocks) / m) ** 2
    se = np.sqrt(np.var(pseudo, axis=0, ddof=1) / denom)

    # Genomic control uses the median chi-square over kept SNPs.
    chi_masked = np.where(keep, chi2, np.nan)
    lambda_gc = np.nanmedian(chi_masked, axis=0) / _CHI2_1DF_MEDIAN

    return BatchedH2Result(
        h2=slope_full / n * m,
        h2_se=se,
        intercept=intercept_full,
        mean_chi2=mean_chi,
        lambda_gc=lambda_gc,
        n_snps=cnt.astype(np.int64),
    )


@frozen
class ExactH2Result:
    h2: float
    h2_se: float
    intercept: float


def exact_h2_single(
    chi2: np.ndarray,
    ld: np.ndarray,
    n: float,
    m: float,
    *,
    n_blocks: int = DEFAULT_N_BLOCKS,
    exclude: np.ndarray | None = None,
) -> ExactH2Result:
    """Reference single-protein heritability via the repo's exact GenomicSEM port. Drops
    (rather than zero-weights) filtered/cis SNPs, then blocks the KEPT SNPs -- this is what
    the batched kernel is validated against. chi2/ld are genome-sorted (S,) arrays."""
    keep = np.isfinite(chi2) & (chi2 <= max(0.001 * n, 80.0))
    if exclude is not None:
        keep &= ~exclude
    chi = chi2[keep]
    ld_kept = ld[keep]
    est = _estimate_h2(
        chi=chi,
        ld_raw=ld_kept,
        wld_raw=ld_kept,
        n=np.full(chi.shape, n),
        m=m,
        n_blocks=n_blocks,
    )
    # run_ldsc: V_h2 = var(pseudo, ddof=1) / (n_bar * sqrt(n_blocks) / m)^2 ; se = sqrt(V).
    denom = (est.n_bar * np.sqrt(n_blocks) / m) ** 2
    se = float(np.sqrt(np.var(est.pseudo_coef, ddof=1) / denom))
    return ExactH2Result(h2=est.reg_tot, h2_se=se, intercept=est.intercept)
