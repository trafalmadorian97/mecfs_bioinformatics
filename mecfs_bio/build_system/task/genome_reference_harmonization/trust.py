"""
Whether a table's NEA column is trustworthy as the reference allele.

Trust requires exactly 100% consistency. That is, every checkable SNV (non-palindromic,
single base) has NEA as the reference base, and every checkable indel (different
lengths) has only NEA matching the reference. Each set must also meet a minimum count.
Ambiguous indels, whose alleles both match, carry no evidence and are not counted.
"""

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.genome_reference_harmonization.allele_classes import (
    ALLELE_CLASS_COL,
    CLASS_INDEL_EA_ONLY,
    CLASS_INDEL_NEA_ONLY,
    CLASS_INDEL_NOT_ON_REFERENCE,
    CLASS_NEA_REF,
    IS_PALINDROMIC_SNV_COL,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
)

_CONSISTENT_SNVS = "consistent_snvs"
_INCONSISTENT_SNVS = "inconsistent_snvs"
_CONSISTENT_INDELS = "consistent_indels"
_INCONSISTENT_INDELS = "inconsistent_indels"


@frozen
class TrustCounts:
    consistent_snvs: int
    inconsistent_snvs: int
    consistent_indels: int
    inconsistent_indels: int

    def __add__(self, other: "TrustCounts") -> "TrustCounts":
        return TrustCounts(
            consistent_snvs=self.consistent_snvs + other.consistent_snvs,
            inconsistent_snvs=self.inconsistent_snvs + other.inconsistent_snvs,
            consistent_indels=self.consistent_indels + other.consistent_indels,
            inconsistent_indels=self.inconsistent_indels + other.inconsistent_indels,
        )

    @classmethod
    def zero(cls) -> "TrustCounts":
        return cls(
            consistent_snvs=0,
            inconsistent_snvs=0,
            consistent_indels=0,
            inconsistent_indels=0,
        )


def count_trust_evidence(classified: pl.DataFrame) -> TrustCounts:
    """Count checkable SNVs and indels in a frame produced by classify_alleles."""
    allele_class = pl.col(ALLELE_CLASS_COL)
    ea_length = pl.col(GWASLAB_EFFECT_ALLELE_COL).str.len_bytes()
    nea_length = pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).str.len_bytes()
    checkable_snv = (
        (ea_length == 1) & (nea_length == 1) & ~pl.col(IS_PALINDROMIC_SNV_COL)
    )
    indel = ea_length != nea_length
    counts = classified.select(
        (checkable_snv & (allele_class == CLASS_NEA_REF)).sum().alias(_CONSISTENT_SNVS),
        (checkable_snv & (allele_class != CLASS_NEA_REF))
        .sum()
        .alias(_INCONSISTENT_SNVS),
        (indel & (allele_class == CLASS_INDEL_NEA_ONLY))
        .sum()
        .alias(_CONSISTENT_INDELS),
        (
            indel
            & allele_class.is_in([CLASS_INDEL_EA_ONLY, CLASS_INDEL_NOT_ON_REFERENCE])
        )
        .sum()
        .alias(_INCONSISTENT_INDELS),
    ).row(0, named=True)
    return TrustCounts(
        consistent_snvs=int(counts[_CONSISTENT_SNVS]),
        inconsistent_snvs=int(counts[_INCONSISTENT_SNVS]),
        consistent_indels=int(counts[_CONSISTENT_INDELS]),
        inconsistent_indels=int(counts[_INCONSISTENT_INDELS]),
    )


def decide_trust(
    counts: TrustCounts, options: GenomeReferenceHarmonizationOptions
) -> bool:
    return (
        counts.inconsistent_snvs == 0
        and counts.inconsistent_indels == 0
        and counts.consistent_snvs >= options.min_checkable_snvs
        and counts.consistent_indels >= options.min_checkable_indels
    )
