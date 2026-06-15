"""
Reference-align and standardise per-SNP GWAS effects for GWAS-by-subtraction.

Each trait's effects are aligned to a reference panel and rescaled to a single
linear (OLS-style) standardised scale, then all traits are listwise-merged onto
the (MAF-filtered) reference. The result is the per-SNP table the kernel
consumes: SNP, CHR, BP, MAF, A1, A2, and beta.<trait>/se.<trait> columns. (MAF
here is the *reference* MAF; the kernel derives varSNP from it.)

**Every trait is treated linearly.** Unlike ``GenomicSEM::sumstats``, there is
no per-trait method flag and no binary/logistic (``se.logit``/``linprob``)
standardisation. GenomicSEM standardises a binary trait's effects to the
logistic-latent scale (``beta = logOR/sqrt(logOR^2*varSNP + pi^2/3)``) while
``ldsc`` puts the genetic-covariance matrix on the prevalence-converted
liability scale; combining those two scales in the subtraction leaves the
remainder factor non-orthogonal to the reference trait. We deliberately drop
that machinery and standardise everything as ``beta = Z/sqrt(N*varSNP)`` so the
genetic-covariance matrix and the per-SNP betas share one scale. See
``genomic_sem_gwas_by_subtraction_full_python_task`` for the full rationale.

Per trait, in order:

1. Override N with the provided scalar (when given).
2. Drop SNPs that appear more than once (multiallelic).
3. Upper-case A1/A2, null out non-ACGT (keep_indel=False).
4. Inner-merge with the reference on SNP (.x = reference, .y = file).
5. Drop rows with missing P or effect.
6. varSNP from the *file* MAF when present (folded to minor, dropping 0/1),
   else from the reference MAF.
7. Odds-ratio guard: if round(median(effect)) == 1 the effect column looks
   like an OR. The standardisation uses sign(effect), so an OR-coded column
   (always > 0) would make every Z positive; raise rather than proceed.
8. Drop effect == 0.
9. Z = sign(effect) * Phi^{-1}(1 - P/2); for P < 1e-307, sign(effect)*sqrt(-2 ln P).
10. Standardise: effect = Z / sqrt(N * varSNP).
11. Flip effect to the reference A1 allele; drop allele mismatches.
12. INFO filter when present.
13. Emit (SNP, beta=effect, se=|effect/Z|).
"""

from __future__ import annotations

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
    n: float | None  # provided total sample size (None -> use the file's N column)


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
    """Return a DataFrame with columns SNP, beta, se for one trait.

    Note:
        - In this method, the ref dataframe is typically 1000 genomes
        - This is in contrast to munging, where the ref dataframe is hapmap3
    """
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

    # Z, with a large-deviation approximation for extremely small P. Below
    # _TINY_P, p/2 underflows in double precision and norm.isf saturates to inf,
    # so use z ~ sqrt(-2 ln p) there instead. P is already finite and non-null
    # here (filtered above), so no NaN/inf guard is needed.
    p = merged[MUNGE_P_COL].to_numpy().astype(float)
    effect = merged[MUNGE_EFFECT_COL].to_numpy().astype(float)
    tiny = p < _TINY_P
    z = np.empty_like(p)
    z[~tiny] = np.sign(effect[~tiny]) * norm.isf(p[~tiny] / 2.0)
    if tiny.any():
        z[tiny] = np.sign(effect[tiny]) * np.sqrt(-2.0 * np.log(p[tiny]))
    merged = merged.with_columns(pl.Series(MUNGE_Z_COL, z))

    # Linear (OLS-style) standardisation of the per-SNP effect: beta on the
    # standardised-genotype, standardised-phenotype scale. Applied to every trait
    # -- there is intentionally no binary/logistic special-casing (see the module
    # docstring). See notes on Sampling noise and LDSC.
    merged = merged.with_columns(
        (
            pl.col(MUNGE_Z_COL) / (pl.col(MUNGE_N_COL) * pl.col(_VARSNP_COL)).sqrt()
        ).alias(MUNGE_EFFECT_COL)
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

    # Output: linear-scale beta and its SE (se = |beta / Z|).
    out = merged.select(
        pl.col(MUNGE_SNP_COL),
        pl.col(MUNGE_EFFECT_COL).alias(_BETA_OUT_COL),
        (pl.col(MUNGE_EFFECT_COL) / pl.col(MUNGE_Z_COL)).abs().alias(_SE_OUT_COL),
    )
    # R's na.omit drops rows with NA *or* NaN. polars drop_nulls keeps NaN, so
    # filter both: e.g. a SNP with P == 1 gives Z == 0 -> effect == 0 ->
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
