"""
Pure-Python (polars) re-implementation of ``GenomicSEM::munge``.

``munge`` QC's a single trait's GWAS summary statistics against a HapMap3
reference and produces the ``.sumstats`` table (SNP, N, Z, A1, A2) that
:func:`_genomic_sem_ldsc.run_ldsc` consumes. This port operates on a polars
DataFrame so the full pipeline can munge straight from the in-memory source
(no multi-minute R read of giant text files).

Input columns (canonical munge names, as written by ``_write_munge_input``):
SNP, A1, A2, effect, P, N, and optionally MAF / INFO. (SE is accepted but
unused, matching R: Z is derived from the effect sign and P.)

The steps mirror ``GenomicSEM:::.munge_main`` exactly, in the same order:

1. Override N with the provided scalar (when given).
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

_ACGT = ["A", "C", "G", "T"]


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
    assert "P" in df.columns, "munge input requires a P column"
    assert "effect" in df.columns, "munge input requires an effect column"
    has_maf = "MAF" in df.columns
    has_info = "INFO" in df.columns

    work = df
    if n is not None:
        work = work.with_columns(pl.lit(float(n)).alias("N"))
    if has_maf:
        work = work.with_columns(
            pl.when(pl.col("MAF") <= 0.5)
            .then(pl.col("MAF"))
            .otherwise(1.0 - pl.col("MAF"))
            .alias("MAF")
        )
    work = work.with_columns(_restrict_to_acgt("A1"), _restrict_to_acgt("A2"))

    ref_aligned = ref.select(
        pl.col("SNP"),
        pl.col("A1").str.to_uppercase().alias("A1_ref"),
        pl.col("A2").str.to_uppercase().alias("A2_ref"),
    )
    merged = ref_aligned.join(work, on="SNP", how="inner")
    merged = merged.filter(pl.col("P").is_not_null() & pl.col("effect").is_not_null())

    # Odds-ratio detection on the merged effect column.
    median_effect = merged.select(pl.col("effect").median()).item()
    if median_effect is not None and round(median_effect) == 1:
        merged = merged.with_columns(pl.col("effect").log().alias("effect"))

    # Flip effect to the reference A1 allele.
    merged = merged.with_columns(
        pl.when((pl.col("A1_ref") != pl.col("A1")) & (pl.col("A1_ref") == pl.col("A2")))
        .then(-pl.col("effect"))
        .otherwise(pl.col("effect"))
        .alias("effect")
    )

    # Keep only rows whose alleles match the reference (either orientation).
    # polars .filter drops nulls, matching R's NA-dropping subset().
    merged = merged.filter(
        ((pl.col("A1_ref") == pl.col("A1")) | (pl.col("A1_ref") == pl.col("A2")))
        & ((pl.col("A2_ref") == pl.col("A2")) | (pl.col("A2_ref") == pl.col("A1")))
    )

    if has_info:
        merged = merged.filter(pl.col("INFO") >= info_filter)
    if has_maf:
        merged = merged.filter(
            (pl.col("MAF") >= maf_filter) & pl.col("MAF").is_not_null()
        )

    effect = merged["effect"].to_numpy()
    p = merged["P"].to_numpy()
    z = np.sign(effect) * norm.isf(p / 2.0)

    return merged.select(
        pl.col("SNP"),
        pl.col("N"),
        pl.Series("Z", z),
        pl.col("A1_ref").alias("A1"),
        pl.col("A2_ref").alias("A2"),
    )
