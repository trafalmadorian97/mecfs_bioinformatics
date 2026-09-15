"""
Resolve one chromosome of summary statistics against the genome reference.

Every input row is returned. Alleles are oriented so that NEA is the plus-strand
reference allele, allele-dependent statistics are flipped to match, and
DROP_REASON_COL is null for rows to keep. The function holds one chromosome in memory.
"""

from collections.abc import Mapping
from typing import Protocol

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.genome_reference_harmonization.allele_classes import (
    ALLELE_CLASS_COL,
    CLASS_EA_REF,
    CLASS_EA_REF_RC,
    CLASS_INDEL_BOTH,
    CLASS_INDEL_EA_ONLY,
    CLASS_INDEL_NOT_ON_REFERENCE,
    CLASS_NEA_REF_RC,
    CLASS_NOT_ON_REFERENCE,
    IS_PALINDROMIC_MNP_COL,
    IS_PALINDROMIC_SNV_COL,
    classify_alleles,
    prepare_alleles,
    reverse_complement_expr,
    valid_alleles_expr,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.flip import (
    FlipRule,
    flip_statistics,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.outcomes import (
    ACTION_COMPLEMENT,
    ACTION_COMPLEMENT_SWAP,
    ACTION_KEEP,
    ACTION_SWAP,
    DROP_INDEL_NOT_ON_REFERENCE,
    DROP_INVALID_ALLELE,
    DROP_NOT_ON_REFERENCE,
    DROP_PALINDROME_UNRESOLVED,
    DROP_PALINDROMIC_MNP_UNTRUSTED,
    PALINDROME_STRAND_FLIP,
    PALINDROME_UNRESOLVED,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.palindromes import (
    decide_palindrome_strands,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

DROP_REASON_COL = "_drop_reason"
ROW_INDEX_COL = "_row_index"
ACTION_COL = "_action"
PALINDROME_DECISION_COL = "_palindrome_decision"


class PanelLoader(Protocol):
    def __call__(self, positions: pl.Series) -> pl.DataFrame:
        """Panel rows (POS Int64, REF, ALT, AF) of the current chromosome at these positions."""
        ...


@frozen
class ChromosomeContext:
    chrom: int
    fasta: IndexedFasta
    trusted: bool
    rules: Mapping[str, FlipRule]
    options: GenomeReferenceHarmonizationOptions


def resolve_chromosome(
    rows: pl.DataFrame, context: ChromosomeContext, load_panel: PanelLoader
) -> pl.DataFrame:
    output_columns = [*rows.columns, DROP_REASON_COL]
    frame = prepare_alleles(rows).with_row_index(ROW_INDEX_COL)
    invalid = frame.filter(~valid_alleles_expr()).with_columns(
        pl.lit(DROP_INVALID_ALLELE, dtype=pl.String).alias(DROP_REASON_COL)
    )
    valid = classify_alleles(
        frame.filter(valid_alleles_expr()),
        fasta=context.fasta,
        chrom=context.chrom,
        max_gather_bytes=context.options.max_gather_bytes,
    ).with_columns(_base_action_expr(), _base_drop_reason_expr(trusted=context.trusted))
    panel = None if context.trusted else load_panel(_panel_positions(valid))
    valid = _apply_allele_actions(valid, context.rules)
    if panel is not None:
        valid = _apply_palindrome_strand_rules(valid, panel, context)
    return (
        pl.concat(
            [
                valid.select(ROW_INDEX_COL, *output_columns),
                invalid.select(ROW_INDEX_COL, *output_columns),
            ]
        )
        .sort(ROW_INDEX_COL)
        .drop(ROW_INDEX_COL)
    )


def _base_action_expr() -> pl.Expr:
    allele_class = pl.col(ALLELE_CLASS_COL)
    return (
        pl.when(allele_class.is_in([CLASS_EA_REF, CLASS_INDEL_EA_ONLY]))
        .then(pl.lit(ACTION_SWAP))
        .when(allele_class == CLASS_NEA_REF_RC)
        .then(pl.lit(ACTION_COMPLEMENT))
        .when(allele_class == CLASS_EA_REF_RC)
        .then(pl.lit(ACTION_COMPLEMENT_SWAP))
        .otherwise(pl.lit(ACTION_KEEP))
        .alias(ACTION_COL)
    )


def _base_drop_reason_expr(trusted: bool) -> pl.Expr:
    allele_class = pl.col(ALLELE_CLASS_COL)
    untrusted_palindromic_mnp = pl.col(IS_PALINDROMIC_MNP_COL) & pl.lit(not trusted)
    return (
        pl.when(allele_class == CLASS_NOT_ON_REFERENCE)
        .then(pl.lit(DROP_NOT_ON_REFERENCE))
        .when(allele_class == CLASS_INDEL_NOT_ON_REFERENCE)
        .then(pl.lit(DROP_INDEL_NOT_ON_REFERENCE))
        .when(untrusted_palindromic_mnp)
        .then(pl.lit(DROP_PALINDROMIC_MNP_UNTRUSTED))
        .otherwise(pl.lit(None, dtype=pl.String))
        .alias(DROP_REASON_COL)
    )


def _apply_allele_actions(
    frame: pl.DataFrame, rules: Mapping[str, FlipRule]
) -> pl.DataFrame:
    ea = pl.col(GWASLAB_EFFECT_ALLELE_COL)
    nea = pl.col(GWASLAB_NON_EFFECT_ALLELE_COL)
    complement = pl.col(ACTION_COL).is_in([ACTION_COMPLEMENT, ACTION_COMPLEMENT_SWAP])
    swap = pl.col(ACTION_COL).is_in([ACTION_SWAP, ACTION_COMPLEMENT_SWAP])
    complemented = frame.with_columns(
        pl.when(complement)
        .then(reverse_complement_expr(GWASLAB_EFFECT_ALLELE_COL))
        .otherwise(ea)
        .alias(GWASLAB_EFFECT_ALLELE_COL),
        pl.when(complement)
        .then(reverse_complement_expr(GWASLAB_NON_EFFECT_ALLELE_COL))
        .otherwise(nea)
        .alias(GWASLAB_NON_EFFECT_ALLELE_COL),
    )
    swapped = complemented.with_columns(
        pl.when(swap).then(nea).otherwise(ea).alias(GWASLAB_EFFECT_ALLELE_COL),
        pl.when(swap).then(ea).otherwise(nea).alias(GWASLAB_NON_EFFECT_ALLELE_COL),
    )
    return flip_statistics(swapped, mask=swap, rules=rules)


def _panel_positions(valid: pl.DataFrame) -> pl.Series:
    """Positions whose variants may need panel evidence: palindromic SNVs and ambiguous indels."""
    needs_panel = pl.col(IS_PALINDROMIC_SNV_COL) | (
        pl.col(ALLELE_CLASS_COL) == CLASS_INDEL_BOTH
    )
    return valid.filter(needs_panel)[GWASLAB_POS_COL]


def _eaf_expr(frame: pl.DataFrame) -> pl.Expr:
    if GWASLAB_EFFECT_ALLELE_FREQ_COL in frame.columns:
        return pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL).cast(pl.Float64)
    return pl.lit(None, dtype=pl.Float64).alias(GWASLAB_EFFECT_ALLELE_FREQ_COL)


def _apply_palindrome_strand_rules(
    valid: pl.DataFrame, panel: pl.DataFrame, context: ChromosomeContext
) -> pl.DataFrame:
    candidates = valid.filter(
        pl.col(IS_PALINDROMIC_SNV_COL) & pl.col(DROP_REASON_COL).is_null()
    )
    decisions = candidates.select(ROW_INDEX_COL).with_columns(
        decide_palindrome_strands(
            candidates.select(
                GWASLAB_POS_COL,
                GWASLAB_EFFECT_ALLELE_COL,
                GWASLAB_NON_EFFECT_ALLELE_COL,
                _eaf_expr(candidates),
            ),
            panel,
            context.options,
        ).alias(PALINDROME_DECISION_COL)
    )
    decided = valid.join(decisions, on=ROW_INDEX_COL, how="left", maintain_order="left")
    decided = flip_statistics(
        decided,
        mask=pl.col(PALINDROME_DECISION_COL) == PALINDROME_STRAND_FLIP,
        rules=context.rules,
    )
    if not context.options.keep_unresolved_palindromes:
        unresolved = (
            pl.col(PALINDROME_DECISION_COL) == PALINDROME_UNRESOLVED
        ) & pl.col(DROP_REASON_COL).is_null()
        decided = decided.with_columns(
            pl.when(unresolved)
            .then(pl.lit(DROP_PALINDROME_UNRESOLVED))
            .otherwise(pl.col(DROP_REASON_COL))
            .alias(DROP_REASON_COL)
        )
    return decided.drop(PALINDROME_DECISION_COL)
