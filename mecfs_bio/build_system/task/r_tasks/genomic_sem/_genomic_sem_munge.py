"""
Pure-Python (polars) re-implementation of ``GenomicSEM::munge``.

``munge`` QC's a single trait's GWAS summary statistics against a HapMap3
reference and produces the ``.sumstats`` table (SNP, N, Z, A1, A2) that
:func:`_genomic_sem_ldsc.run_ldsc` consumes. This port operates on a polars
DataFrame so the full pipeline can munge straight from the in-memory source
(no multi-minute R read of giant text files).

Input columns (canonical munge names, as written by ``write_munge_input``):
SNP, A1, A2, effect, P, N, and optionally MAF / INFO. (SE is accepted but
unused, matching R: Z is derived from the effect sign and P.)

The steps mirror ``GenomicSEM:::.munge_main`` exactly, in the same order:

1. Override N with the provided scalar (when given and not NaN; a NaN N means
   "not provided", matching R's `!is.na(N)`, so the file's N column is kept).
2. Fold MAF to the minor-allele frequency (min(MAF, 1-MAF)).
3. Upper-case A1/A2 and null out anything outside {A, C, G, T}.
4. Inner-merge with the reference on SNP (reference alleles become A1.x/A2.x).
5. Drop rows with missing P or effect.
6. Odds-ratio detection: if round(median(effect)) == 1, take log(effect).
7. Flip the effect sign where the file's effect allele matches the
   reference's *other* allele.
8. Drop rows whose alleles don't match the reference (either orientation).
9. Z = sign(effect) * sqrt(qchisq(P, 1, lower=FALSE))
        = sign(effect) * Phi^{-1}(1 - P/2).
10. Filter on INFO >= info_filter and MAF >= maf_filter when present.

Output: SNP, N, Z, A1, A2 where A1/A2 are the *reference* alleles (so Z is the
effect of the reference A1 allele).
"""

from __future__ import annotations

import numpy as np
import polars as pl
from scipy.stats import norm

from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_config import (
    MUNGE_A1_COL,
    MUNGE_A2_COL,
    MUNGE_EFFECT_COL,
    MUNGE_INFO_COL,
    MUNGE_MAF_COL,
    MUNGE_N_COL,
    MUNGE_P_COL,
    MUNGE_SNP_COL,
    MUNGE_Z_COL,
)

_ACGT = ["A", "C", "G", "T"]

# Internal join-suffix names for the reference panel's alleles (not canonical
# munge columns), so the file's A1/A2 can be compared against them.
_A1_REF_COL = f"{MUNGE_A1_COL}_ref"
_A2_REF_COL = f"{MUNGE_A2_COL}_ref"


def _restrict_to_acgt(col: str) -> pl.Expr:
    upper = pl.col(col).str.to_uppercase()
    return pl.when(upper.is_in(_ACGT)).then(upper).otherwise(None).alias(col)


def munge_sumstats(
    df: pl.DataFrame,
    ref: pl.DataFrame,
    *,
    n: float | None,
    info_filter: float = 0.9,
    maf_filter: float = 0.01,
) -> pl.DataFrame:
    """
    Munge one trait. ``df`` holds canonical munge columns; ``ref`` is the
    HapMap3 snplist (SNP, A1, A2). Returns a polars DataFrame with columns
    SNP, N, Z, A1, A2.
    """
    assert MUNGE_P_COL in df.columns, "munge input requires a P column"
    assert MUNGE_EFFECT_COL in df.columns, "munge input requires an effect column"
    has_maf = MUNGE_MAF_COL in df.columns
    has_info = MUNGE_INFO_COL in df.columns

    work = df
    # Override N with the provided scalar only when it is a real number. This
    # mirrors GenomicSEM::munge's `!is.na(N)` guard: a NaN sample size means
    # "not provided", so the file's own N column is kept rather than clobbered.
    if n is not None and not np.isnan(n):
        work = work.with_columns(pl.lit(float(n)).alias(MUNGE_N_COL))
    if has_maf:
        work = work.with_columns(
            pl.when(pl.col(MUNGE_MAF_COL) <= 0.5)
            .then(pl.col(MUNGE_MAF_COL))
            .otherwise(1.0 - pl.col(MUNGE_MAF_COL))
            .alias(MUNGE_MAF_COL)
        )
    work = work.with_columns(
        _restrict_to_acgt(MUNGE_A1_COL), _restrict_to_acgt(MUNGE_A2_COL)
    )

    ref_aligned = ref.select(
        pl.col(MUNGE_SNP_COL),
        pl.col(MUNGE_A1_COL).str.to_uppercase().alias(_A1_REF_COL),
        pl.col(MUNGE_A2_COL).str.to_uppercase().alias(_A2_REF_COL),
    )
    merged = ref_aligned.join(work, on=MUNGE_SNP_COL, how="inner")
    merged = merged.filter(
        pl.col(MUNGE_P_COL).is_not_null() & pl.col(MUNGE_EFFECT_COL).is_not_null()
    )

    # Odds-ratio detection on the merged effect column.
    median_effect = merged.select(pl.col(MUNGE_EFFECT_COL).median()).item()
    if median_effect is not None and round(median_effect) == 1:
        merged = merged.with_columns(
            pl.col(MUNGE_EFFECT_COL).log().alias(MUNGE_EFFECT_COL)
        )

    # Flip effect to the reference A1 allele.
    merged = merged.with_columns(
        pl.when(
            (pl.col(_A1_REF_COL) != pl.col(MUNGE_A1_COL))
            & (pl.col(_A1_REF_COL) == pl.col(MUNGE_A2_COL))
        )
        .then(-pl.col(MUNGE_EFFECT_COL))
        .otherwise(pl.col(MUNGE_EFFECT_COL))
        .alias(MUNGE_EFFECT_COL)
    )

    # Keep only rows whose alleles match the reference (either orientation).
    # polars .filter drops nulls, matching R's NA-dropping subset().
    merged = merged.filter(
        (
            (pl.col(_A1_REF_COL) == pl.col(MUNGE_A1_COL))
            | (pl.col(_A1_REF_COL) == pl.col(MUNGE_A2_COL))
        )
        & (
            (pl.col(_A2_REF_COL) == pl.col(MUNGE_A2_COL))
            | (pl.col(_A2_REF_COL) == pl.col(MUNGE_A1_COL))
        )
    )

    if has_info:
        merged = merged.filter(pl.col(MUNGE_INFO_COL) >= info_filter)
    if has_maf:
        merged = merged.filter(
            (pl.col(MUNGE_MAF_COL) >= maf_filter) & pl.col(MUNGE_MAF_COL).is_not_null()
        )

    effect = merged[MUNGE_EFFECT_COL].to_numpy()
    p = merged[MUNGE_P_COL].to_numpy()
    z = np.sign(effect) * norm.isf(p / 2.0)

    return merged.select(
        pl.col(MUNGE_SNP_COL),
        pl.col(MUNGE_N_COL),
        pl.Series(MUNGE_Z_COL, z),
        pl.col(_A1_REF_COL).alias(MUNGE_A1_COL),
        pl.col(_A2_REF_COL).alias(MUNGE_A2_COL),
    )
