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
7. Odds-ratio guard: if round(median(effect)) == 1 the effect column looks
   like an OR, which this pipeline never produces (gwaslab keeps ORs separate
   and feeds only BETA into `effect`), so raise rather than silently log() as
   GenomicSEM does.
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

from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_config import (
    MUNGE_A1_COL,
    MUNGE_A2_COL,
    MUNGE_EFFECT_COL,
    MUNGE_INFO_COL,
    MUNGE_MAF_COL,
    MUNGE_N_COL,
    MUNGE_P_COL,
    MUNGE_SE_COL,
    MUNGE_SNP_COL,
    MUNGE_Z_COL,
)

_ACGT = ["A", "C", "G", "T"]
_TINY_P = 1e-307

# Internal column names introduced while aligning a trait to the reference
# (the reference side gets a "_ref" suffix, the file side a "_file" suffix);
# these are not canonical munge columns. varSNP is the derived per-SNP variance.
_A1_REF_COL = f"{MUNGE_A1_COL}_ref"
_A2_REF_COL = f"{MUNGE_A2_COL}_ref"
_A1_FILE_COL = f"{MUNGE_A1_COL}_file"
_A2_FILE_COL = f"{MUNGE_A2_COL}_file"
_MAF_REF_COL = "maf_ref"
_MAF_FILE_COL = "maf_file"
_VARSNP_COL = "varSNP"

# Per-trait output columns; run_sumstats prefixes these with the trait name
# (e.g. "beta.<trait>" / "se.<trait>").
_BETA_OUT_COL = "beta"
_SE_OUT_COL = "se"


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
    ref: pl.DataFrame, maf_filter: float, exclude_ambig: bool
) -> pl.DataFrame:
    out = ref.filter(pl.col(MUNGE_MAF_COL) >= maf_filter)
    if exclude_ambig:
        a1 = pl.col(MUNGE_A1_COL).str.to_uppercase()
        a2 = pl.col(MUNGE_A2_COL).str.to_uppercase()
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
    assert MUNGE_P_COL in work.columns and MUNGE_EFFECT_COL in work.columns
    # Override N with the provided scalar only when it is a real number (mirrors
    # GenomicSEM's `!is.na(N)`): a NaN sample size means "not provided", so the
    # file's own N column is kept rather than clobbered.
    if trait.n is not None and not np.isnan(trait.n):
        work = work.with_columns(pl.lit(float(trait.n)).alias(MUNGE_N_COL))

    # Drop multiallelic (duplicated) SNPs entirely.
    work = work.filter(pl.col(MUNGE_SNP_COL).is_unique())
    work = work.with_columns(
        _restrict_to_acgt(MUNGE_A1_COL), _restrict_to_acgt(MUNGE_A2_COL)
    )

    ref_aligned = ref.rename(
        {
            MUNGE_A1_COL: _A1_REF_COL,
            MUNGE_A2_COL: _A2_REF_COL,
            MUNGE_MAF_COL: _MAF_REF_COL,
        }
    )
    has_file_maf = MUNGE_MAF_COL in work.columns
    if has_file_maf:
        work = work.rename({MUNGE_MAF_COL: _MAF_FILE_COL})
    work = work.rename({MUNGE_A1_COL: _A1_FILE_COL, MUNGE_A2_COL: _A2_FILE_COL})

    merged = ref_aligned.join(work, on=MUNGE_SNP_COL, how="inner")
    # Drop rows with missing P or effect. R's is.na() is true for NA *and* NaN,
    # but polars is_not_null keeps NaN, so we also exclude NaN to match R (and to
    # keep the odds-ratio median below from being poisoned by a stray NaN).
    merged = merged.filter(
        pl.col(MUNGE_P_COL).is_not_null()
        & pl.col(MUNGE_P_COL).is_not_nan()
        & pl.col(MUNGE_EFFECT_COL).is_not_null()
        & pl.col(MUNGE_EFFECT_COL).is_not_nan()
    )

    # varSNP from the file MAF (folded, dropping monomorphic) or the ref MAF.
    if has_file_maf:
        merged = merged.with_columns(
            pl.when(pl.col(_MAF_FILE_COL) > 0.5)
            .then(1.0 - pl.col(_MAF_FILE_COL))
            .otherwise(pl.col(_MAF_FILE_COL))
            .alias(_MAF_FILE_COL)
        )
        merged = merged.filter(
            (pl.col(_MAF_FILE_COL) != 0) & (pl.col(_MAF_FILE_COL) != 1)
        )
        merged = merged.with_columns(
            (2.0 * pl.col(_MAF_FILE_COL) * (1.0 - pl.col(_MAF_FILE_COL))).alias(
                _VARSNP_COL
            )
        )
    else:
        merged = merged.with_columns(
            (2.0 * pl.col(_MAF_REF_COL) * (1.0 - pl.col(_MAF_REF_COL))).alias(
                _VARSNP_COL
            )
        )

    # Guard against an odds-ratio effect column. GenomicSEM silently
    # log-transforms when round(median(effect)) == 1 (assuming the column is an
    # OR), but this pipeline always feeds a log-scale beta -- gwaslab keeps odds
    # ratios in a separate OR column, and build_munge_input_df wires only BETA
    # into `effect`. A median near 1 therefore signals a misspecified input
    # rather than an OR to convert, so fail loudly instead of mutating the data.
    median_effect = merged.select(pl.col(MUNGE_EFFECT_COL).median()).item()
    if median_effect is not None and round(median_effect) == 1:
        raise ValueError(
            f"effect column has median {median_effect:.4g} (rounds to 1), which "
            "looks like an odds ratio. This pipeline expects log-scale betas; "
            "supply odds ratios via the gwaslab OR column / convert to beta "
            "upstream rather than passing them as 'effect'."
        )
    merged = merged.filter(pl.col(MUNGE_EFFECT_COL) != 0)

    if merged.height == 0:
        return pl.DataFrame({MUNGE_SNP_COL: [], _BETA_OUT_COL: [], _SE_OUT_COL: []})

    # Z with a high-magnitude approximation for extremely small P.
    p = merged[MUNGE_P_COL].to_numpy().astype(float)
    effect = merged[MUNGE_EFFECT_COL].to_numpy().astype(float)
    tiny = ~np.isfinite(p) | (p < _TINY_P)
    z = np.empty_like(p)
    z[~tiny] = np.sign(effect[~tiny]) * norm.isf(p[~tiny] / 2.0)
    if tiny.any():
        z[tiny] = np.sign(effect[tiny]) * np.sqrt(-2.0 * np.log(p[tiny]))
    merged = merged.with_columns(pl.Series(MUNGE_Z_COL, z))

    # Method-specific rescaling of effect (and SE for linprob).
    if trait.ols:
        merged = merged.with_columns(
            (
                pl.col(MUNGE_Z_COL) / (pl.col(MUNGE_N_COL) * pl.col(_VARSNP_COL)).sqrt()
            ).alias(MUNGE_EFFECT_COL)
        )
    elif trait.linprob:
        merged = merged.with_columns(
            (
                pl.col(MUNGE_Z_COL)
                / ((pl.col(MUNGE_N_COL) / 4.0) * pl.col(_VARSNP_COL)).sqrt()
            ).alias(MUNGE_EFFECT_COL),
            (1.0 / ((pl.col(MUNGE_N_COL) / 4.0) * pl.col(_VARSNP_COL)).sqrt()).alias(
                MUNGE_SE_COL
            ),
        )

    # Flip effect to the reference A1 allele, then drop allele mismatches.
    merged = merged.with_columns(
        pl.when(
            (pl.col(_A1_REF_COL) != pl.col(_A1_FILE_COL))
            & (pl.col(_A1_REF_COL) == pl.col(_A2_FILE_COL))
        )
        .then(-pl.col(MUNGE_EFFECT_COL))
        .otherwise(pl.col(MUNGE_EFFECT_COL))
        .alias(MUNGE_EFFECT_COL)
    )
    merged = merged.filter(
        (
            (pl.col(_A1_REF_COL) == pl.col(_A1_FILE_COL))
            | (pl.col(_A1_REF_COL) == pl.col(_A2_FILE_COL))
        )
        & (
            (pl.col(_A2_REF_COL) == pl.col(_A2_FILE_COL))
            | (pl.col(_A2_REF_COL) == pl.col(_A1_FILE_COL))
        )
    )
    if MUNGE_INFO_COL in merged.columns:
        merged = merged.filter(pl.col(MUNGE_INFO_COL) >= info_filter)

    # Output transform.
    pi_term = (math.pi**2) / 3.0
    if trait.ols:
        out = merged.select(
            pl.col(MUNGE_SNP_COL),
            pl.col(MUNGE_EFFECT_COL).alias(_BETA_OUT_COL),
            (pl.col(MUNGE_EFFECT_COL) / pl.col(MUNGE_Z_COL)).abs().alias(_SE_OUT_COL),
        )
    elif trait.linprob:
        den = (pl.col(MUNGE_EFFECT_COL) ** 2 * pl.col(_VARSNP_COL) + pi_term).sqrt()
        out = merged.select(
            pl.col(MUNGE_SNP_COL),
            (pl.col(MUNGE_EFFECT_COL) / den).alias(_BETA_OUT_COL),
            (pl.col(MUNGE_SE_COL) / den).alias(_SE_OUT_COL),
        ).filter((pl.col(_BETA_OUT_COL) != 0) & (pl.col(_SE_OUT_COL) != 0))
    elif trait.se_logit:
        den = (pl.col(MUNGE_EFFECT_COL) ** 2 * pl.col(_VARSNP_COL) + pi_term).sqrt()
        out = merged.select(
            pl.col(MUNGE_SNP_COL),
            (pl.col(MUNGE_EFFECT_COL) / den).alias(_BETA_OUT_COL),
            (pl.col(MUNGE_SE_COL) / den).alias(_SE_OUT_COL),
        )
    else:
        den = (pl.col(MUNGE_EFFECT_COL) ** 2 * pl.col(_VARSNP_COL) + pi_term).sqrt()
        out = merged.select(
            pl.col(MUNGE_SNP_COL),
            (pl.col(MUNGE_EFFECT_COL) / den).alias(_BETA_OUT_COL),
            (pl.col(MUNGE_SE_COL) / pl.col(MUNGE_EFFECT_COL).exp() / den).alias(
                _SE_OUT_COL
            ),
        )
    # R's na.omit drops rows with NA *or* NaN. polars drop_nulls keeps NaN, so
    # filter both: e.g. an OLS SNP with P == 1 gives Z == 0 -> effect == 0 ->
    # se == |0/0| == NaN, which R discards.
    return out.filter(
        pl.col(_BETA_OUT_COL).is_not_null()
        & pl.col(_SE_OUT_COL).is_not_null()
        & pl.col(_BETA_OUT_COL).is_not_nan()
        & pl.col(_SE_OUT_COL).is_not_nan()
    )


def run_sumstats(
    traits: Sequence[SumstatsTrait],
    ref: pl.DataFrame,
    *,
    maf_filter: float = 0.01,
    info_filter: float = 0.6,
    exclude_ambig: bool = False,
) -> pl.DataFrame:
    """
    Align and standardise all traits against the reference, then listwise-merge.

    Returns a wide DataFrame: the reference columns (SNP, CHR, BP, MAF, A1, A2)
    plus beta.<name> / se.<name> for each trait, restricted to the common SNP
    set. Trait order is preserved (it must match the ldsc trait order).

    ``exclude_ambig=True`` drops strand-ambiguous (A/T, C/G) SNPs from the
    reference before alignment.
    """
    ref_f = _filter_reference(ref, maf_filter, exclude_ambig)
    out = ref_f
    for trait in traits:
        per = _standardize_trait(trait, ref_f, info_filter=info_filter)
        per = per.rename(
            {
                _BETA_OUT_COL: f"{_BETA_OUT_COL}.{trait.name}",
                _SE_OUT_COL: f"{_SE_OUT_COL}.{trait.name}",
            }
        )
        out = out.join(per, on=MUNGE_SNP_COL, how="inner")
    return out
