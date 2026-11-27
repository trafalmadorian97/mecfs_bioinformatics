"""
The effect allele frequency of deCODE rows, from variant_info_iceland.

deCODE's per-protein files give only ImpMAF, the minor allele frequency, which does not say
whether the effect allele is the minor one. variant_info_iceland gives ImpFreqA1 = f, the
frequency of allele A1, with one row per (position, allele); multi-allelic positions have several
rows with the same A1.

Rule, per row with effect allele EA and minor allele frequency m:
  candidates = variant_info rows at the row's (CHR, POS) with A1 == EA;
  a candidate agrees when |min(f, 1 - f) - m| <= EAF_TOLERANCE;
  EAF = f when at least one candidate agrees and the agreeing f values span at most
  EAF_TOLERANCE; otherwise EAF is null (undetermined). EAF is never m itself.
"""

import polars as pl

from mecfs_bio.constants.decode_somascan_constants import (
    DECODE_IMP_MAF_COL,
    EAF_TOLERANCE,
    VI_A1_COL,
    VI_CHR_COL,
    VI_FREQ_COL,
    VI_POS_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_POS_COL,
)

_ROW_COL = "__row__"
_AGREES_COL = "__agrees__"


def resolve_effect_allele_frequency(
    rows: pl.DataFrame, variant_info: pl.DataFrame
) -> pl.Series:
    """
    rows: CHR (chromosome code), POS, EA and ImpMAF.
    variant_info: deduplicated CHR, POS, A1 and ImpFreqA1.
    Returns Float64 EAF aligned to rows, null where the rule above does not determine it.
    """
    f = pl.col(VI_FREQ_COL)
    agreeing_f = f.filter(pl.col(_AGREES_COL))
    keyed = rows.select(
        pl.col(GWASLAB_CHROM_COL).cast(pl.Int8),
        pl.col(GWASLAB_POS_COL).cast(pl.Int64),
        pl.col(GWASLAB_EFFECT_ALLELE_COL),
        pl.col(DECODE_IMP_MAF_COL).cast(pl.Float64),
    ).with_row_index(_ROW_COL)
    candidates = variant_info.select(
        pl.col(VI_CHR_COL).cast(pl.Int8).alias(GWASLAB_CHROM_COL),
        pl.col(VI_POS_COL).cast(pl.Int64).alias(GWASLAB_POS_COL),
        pl.col(VI_A1_COL).alias(GWASLAB_EFFECT_ALLELE_COL),
        f.cast(pl.Float64),
    )
    resolved = (
        keyed.join(
            candidates,
            on=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL],
            how="left",
        )
        .with_columns(
            (
                (pl.min_horizontal(f, 1 - f) - pl.col(DECODE_IMP_MAF_COL)).abs()
                <= EAF_TOLERANCE
            )
            .fill_null(False)
            .alias(_AGREES_COL)
        )
        .group_by(_ROW_COL)
        .agg(
            pl.when(agreeing_f.max() - agreeing_f.min() <= EAF_TOLERANCE)
            .then(agreeing_f.mean())
            .otherwise(None)
            .alias(GWASLAB_EFFECT_ALLELE_FREQ_COL)
        )
        .sort(_ROW_COL)
    )
    assert resolved.height == rows.height
    return resolved[GWASLAB_EFFECT_ALLELE_FREQ_COL].cast(pl.Float64)
