"""
Pure-Python (polars) re-implementation of ``GenomicSEM::sumstats``.

``sumstats`` aligns each trait's per-SNP effects to a reference panel and
rescales them to a unit-variance-phenotype scale, then listwise-merges all
traits onto the (MAF-filtered) reference. The result is the per-SNP table the
GWAS kernels consume: SNP, CHR, BP, MAF, A1, A2, and beta.<trait>/se.<trait>
columns. (MAF here is the *reference* MAF; the kernels derive varSNP from it.)

Port of ``GenomicSEM:::.sumstats_main`` per trait, in order:

1. Override N with the provided scalar (when given).
2. Drop SNPs that appear more than once (multiallelic).
3. Upper-case A1/A2, null out non-ACGT (keep_indel=False).
4. Inner-merge with the reference on SNP (.x = reference, .y = file).
5. Drop rows with missing P or effect.
6. varSNP from the *file* MAF when present (folded to minor, dropping 0/1),
   else from the reference MAF.
7. Odds-ratio detection: if round(median(effect)) == 1, take log(effect).
8. Drop effect == 0.
9. Z = sign(effect) * Phi^{-1}(1 - P/2); for P < 1e-307, sign(effect)*sqrt(-2 ln P).
10. Method-specific rescaling of effect / SE (OLS, linprob, se.logit, default).
11. Flip effect to the reference A1 allele; drop allele mismatches.
12. INFO filter when present.
13. Emit (SNP, beta, se) under the method's transform.

Exactly one of OLS / se_logit / linprob is typically set; if none, the
"default" OR-with-OR-SE transform is used (matching R).
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
import polars as pl
from attrs import frozen
from scipy.stats import norm

_ACGT = ["A", "C", "G", "T"]
_TINY_P = 1e-307


@frozen
class SumstatsTrait:
    """One trait's inputs for :func:`run_sumstats`."""

    df: pl.DataFrame  # canonical cols: SNP, A1, A2, effect, SE, P, N, MAF?, INFO?
    name: str
    n: float | None  # provided sample size (total for OLS, sum-Neff for binary)
    ols: bool = False
    se_logit: bool = False
    linprob: bool = False


def _restrict_to_acgt(col: str) -> pl.Expr:
    upper = pl.col(col).str.to_uppercase()
    return pl.when(upper.is_in(_ACGT)).then(upper).otherwise(None).alias(col)


def _filter_reference(
    ref: pl.DataFrame, maf_filter: float, ambig: bool
) -> pl.DataFrame:
    out = ref.filter(pl.col("MAF") >= maf_filter)
    if ambig:
        a1 = pl.col("A1").str.to_uppercase()
        a2 = pl.col("A2").str.to_uppercase()
        ambiguous = (
            ((a1 == "T") & (a2 == "A"))
            | ((a1 == "A") & (a2 == "T"))
            | ((a1 == "C") & (a2 == "G"))
            | ((a1 == "G") & (a2 == "C"))
        )
        out = out.filter(~ambiguous)
    return out


def _standardize_trait(
    trait: SumstatsTrait, ref: pl.DataFrame, *, info_filter: float
) -> pl.DataFrame:
    """Return a DataFrame with columns SNP, beta, se for one trait."""
    work = trait.df
    assert "P" in work.columns and "effect" in work.columns
    if trait.n is not None:
        work = work.with_columns(pl.lit(float(trait.n)).alias("N"))

    # Drop multiallelic (duplicated) SNPs entirely.
    work = work.filter(pl.col("SNP").is_unique())
    work = work.with_columns(_restrict_to_acgt("A1"), _restrict_to_acgt("A2"))

    ref_aligned = ref.rename({"A1": "A1_ref", "A2": "A2_ref", "MAF": "maf_ref"})
    has_file_maf = "MAF" in work.columns
    if has_file_maf:
        work = work.rename({"MAF": "maf_file"})
    work = work.rename({"A1": "A1_file", "A2": "A2_file"})

    merged = ref_aligned.join(work, on="SNP", how="inner")
    merged = merged.filter(pl.col("P").is_not_null() & pl.col("effect").is_not_null())

    # varSNP from the file MAF (folded, dropping monomorphic) or the ref MAF.
    if has_file_maf:
        merged = merged.with_columns(
            pl.when(pl.col("maf_file") > 0.5)
            .then(1.0 - pl.col("maf_file"))
            .otherwise(pl.col("maf_file"))
            .alias("maf_file")
        )
        merged = merged.filter((pl.col("maf_file") != 0) & (pl.col("maf_file") != 1))
        merged = merged.with_columns(
            (2.0 * pl.col("maf_file") * (1.0 - pl.col("maf_file"))).alias("varSNP")
        )
    else:
        merged = merged.with_columns(
            (2.0 * pl.col("maf_ref") * (1.0 - pl.col("maf_ref"))).alias("varSNP")
        )

    # Odds-ratio detection on the merged effect column.
    median_effect = merged.select(pl.col("effect").median()).item()
    if median_effect is not None and round(median_effect) == 1:
        merged = merged.with_columns(pl.col("effect").log().alias("effect"))
    merged = merged.filter(pl.col("effect") != 0)

    if merged.height == 0:
        return pl.DataFrame({"SNP": [], "beta": [], "se": []})

    # Z with a high-magnitude approximation for extremely small P.
    p = merged["P"].to_numpy().astype(float)
    effect = merged["effect"].to_numpy().astype(float)
    tiny = ~np.isfinite(p) | (p < _TINY_P)
    z = np.empty_like(p)
    z[~tiny] = np.sign(effect[~tiny]) * norm.isf(p[~tiny] / 2.0)
    if tiny.any():
        z[tiny] = np.sign(effect[tiny]) * np.sqrt(-2.0 * np.log(p[tiny]))
    merged = merged.with_columns(pl.Series("Z", z))

    # Method-specific rescaling of effect (and SE for linprob).
    if trait.ols:
        merged = merged.with_columns(
            (pl.col("Z") / (pl.col("N") * pl.col("varSNP")).sqrt()).alias("effect")
        )
    elif trait.linprob:
        merged = merged.with_columns(
            (pl.col("Z") / ((pl.col("N") / 4.0) * pl.col("varSNP")).sqrt()).alias(
                "effect"
            ),
            (1.0 / ((pl.col("N") / 4.0) * pl.col("varSNP")).sqrt()).alias("SE"),
        )

    # Flip effect to the reference A1 allele, then drop allele mismatches.
    merged = merged.with_columns(
        pl.when(
            (pl.col("A1_ref") != pl.col("A1_file"))
            & (pl.col("A1_ref") == pl.col("A2_file"))
        )
        .then(-pl.col("effect"))
        .otherwise(pl.col("effect"))
        .alias("effect")
    )
    merged = merged.filter(
        (
            (pl.col("A1_ref") == pl.col("A1_file"))
            | (pl.col("A1_ref") == pl.col("A2_file"))
        )
        & (
            (pl.col("A2_ref") == pl.col("A2_file"))
            | (pl.col("A2_ref") == pl.col("A1_file"))
        )
    )
    if "INFO" in merged.columns:
        merged = merged.filter(pl.col("INFO") >= info_filter)

    # Output transform.
    pi_term = (math.pi**2) / 3.0
    if trait.ols:
        out = merged.select(
            pl.col("SNP"),
            pl.col("effect").alias("beta"),
            (pl.col("effect") / pl.col("Z")).abs().alias("se"),
        )
    elif trait.linprob:
        den = (pl.col("effect") ** 2 * pl.col("varSNP") + pi_term).sqrt()
        out = merged.select(
            pl.col("SNP"),
            (pl.col("effect") / den).alias("beta"),
            (pl.col("SE") / den).alias("se"),
        ).filter((pl.col("beta") != 0) & (pl.col("se") != 0))
    elif trait.se_logit:
        den = (pl.col("effect") ** 2 * pl.col("varSNP") + pi_term).sqrt()
        out = merged.select(
            pl.col("SNP"),
            (pl.col("effect") / den).alias("beta"),
            (pl.col("SE") / den).alias("se"),
        )
    else:
        den = (pl.col("effect") ** 2 * pl.col("varSNP") + pi_term).sqrt()
        out = merged.select(
            pl.col("SNP"),
            (pl.col("effect") / den).alias("beta"),
            (pl.col("SE") / pl.col("effect").exp() / den).alias("se"),
        )
    # R's na.omit drops rows with NA *or* NaN. polars drop_nulls keeps NaN, so
    # filter both: e.g. an OLS SNP with P == 1 gives Z == 0 -> effect == 0 ->
    # se == |0/0| == NaN, which R discards.
    return out.filter(
        pl.col("beta").is_not_null()
        & pl.col("se").is_not_null()
        & pl.col("beta").is_not_nan()
        & pl.col("se").is_not_nan()
    )


def run_sumstats(
    traits: Sequence[SumstatsTrait],
    ref: pl.DataFrame,
    *,
    maf_filter: float = 0.01,
    info_filter: float = 0.6,
    ambig: bool = False,
) -> pl.DataFrame:
    """
    Align and standardise all traits against the reference, then listwise-merge.

    Returns a wide DataFrame: the reference columns (SNP, CHR, BP, MAF, A1, A2)
    plus beta.<name> / se.<name> for each trait, restricted to the common SNP
    set. Trait order is preserved (it must match the ldsc trait order).
    """
    ref_f = _filter_reference(ref, maf_filter, ambig)
    out = ref_f
    for trait in traits:
        per = _standardize_trait(trait, ref_f, info_filter=info_filter)
        per = per.rename({"beta": f"beta.{trait.name}", "se": f"se.{trait.name}"})
        out = out.join(per, on="SNP", how="inner")
    return out
