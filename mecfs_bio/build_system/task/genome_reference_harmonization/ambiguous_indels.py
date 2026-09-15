"""
Stringent resolution of ambiguous indels (both alleles on the genome) in untrusted tables.

Two readings of a row (EA, NEA, EAF) are possible.

- keep: the variant is REF=NEA, ALT=EA. Its panel record predicts EAF = AF.
- flip: the variant is REF=EA, ALT=NEA. Its panel record predicts EAF = 1 - AF, since EA
  is then the reference allele.

A reading is chosen only if its distance |EAF - predicted| is at most
indel_max_af_distance and, when the other reading also has a record, beats it by at
least indel_min_af_margin. Anything else is dropped: excluding correct variants is
preferred to keeping one whose orientation is wrong.
"""

import polars as pl

from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.outcomes import (
    ACTION_KEEP,
    ACTION_SWAP,
    DROP_AMBIGUOUS_INDEL_AF_INDECISIVE,
    DROP_AMBIGUOUS_INDEL_AF_MISMATCH,
    DROP_AMBIGUOUS_INDEL_NO_EAF,
    DROP_AMBIGUOUS_INDEL_NOT_IN_PANEL,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

INDEL_ACTION_COL = "_indel_action"
INDEL_DROP_REASON_COL = "_indel_drop_reason"
_AF_KEEP = "_af_keep"
_AF_FLIP = "_af_flip"
_KEYS = [GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL]


def _records(
    panel: pl.DataFrame, ref_as: str, alt_as: str, af_name: str
) -> pl.DataFrame:
    return panel.select(
        GWASLAB_POS_COL,
        pl.col(PANEL_REF_COL).alias(ref_as),
        pl.col(PANEL_ALT_COL).alias(alt_as),
        pl.col(PANEL_AF_COL).cast(pl.Float64).alias(af_name),
    )


def decide_ambiguous_indels(
    rows: pl.DataFrame,
    panel: pl.DataFrame,
    options: GenomeReferenceHarmonizationOptions,
) -> pl.DataFrame:
    """INDEL_ACTION_COL and INDEL_DROP_REASON_COL for each row of rows (POS, EA, NEA, EAF), in order."""
    keep_records = _records(
        panel,
        ref_as=GWASLAB_NON_EFFECT_ALLELE_COL,
        alt_as=GWASLAB_EFFECT_ALLELE_COL,
        af_name=_AF_KEEP,
    )
    flip_records = _records(
        panel,
        ref_as=GWASLAB_EFFECT_ALLELE_COL,
        alt_as=GWASLAB_NON_EFFECT_ALLELE_COL,
        af_name=_AF_FLIP,
    )
    joined = rows.join(keep_records, on=_KEYS, how="left", maintain_order="left").join(
        flip_records, on=_KEYS, how="left", maintain_order="left"
    )
    assert joined.height == rows.height, (
        "duplicate panel records for an ambiguous indel"
    )
    eaf = pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL).cast(pl.Float64)
    has_keep = pl.col(_AF_KEEP).is_not_null()
    has_flip = pl.col(_AF_FLIP).is_not_null()
    keep_distance = (eaf - pl.col(_AF_KEEP)).abs()
    flip_distance = (eaf - (1 - pl.col(_AF_FLIP))).abs()
    tolerance = options.indel_max_af_distance
    margin = options.indel_min_af_margin
    choose_keep = (
        has_keep
        & (keep_distance <= tolerance)
        & (~has_flip | (flip_distance - keep_distance >= margin))
    ).fill_null(False)
    choose_flip = (
        has_flip
        & (flip_distance <= tolerance)
        & (~has_keep | (keep_distance - flip_distance >= margin))
    ).fill_null(False)
    return joined.select(
        pl.when(choose_keep)
        .then(pl.lit(ACTION_KEEP))
        .when(choose_flip)
        .then(pl.lit(ACTION_SWAP))
        .otherwise(pl.lit(None, dtype=pl.String))
        .alias(INDEL_ACTION_COL),
        pl.when(eaf.is_null())
        .then(pl.lit(DROP_AMBIGUOUS_INDEL_NO_EAF))
        .when(~has_keep & ~has_flip)
        .then(pl.lit(DROP_AMBIGUOUS_INDEL_NOT_IN_PANEL))
        .when(choose_keep | choose_flip)
        .then(pl.lit(None, dtype=pl.String))
        .when(has_keep ^ has_flip)
        .then(pl.lit(DROP_AMBIGUOUS_INDEL_AF_MISMATCH))
        .otherwise(pl.lit(DROP_AMBIGUOUS_INDEL_AF_INDECISIVE))
        .alias(INDEL_DROP_REASON_COL),
    )
