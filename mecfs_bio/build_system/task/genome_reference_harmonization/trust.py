"""
Whether a table's NEA column is trustworthy as the reference allele.

Trust requires exactly 100% consistency. That is, every checkable SNV (non-palindromic,
single base) has NEA as the reference base, and every checkable indel (different
lengths) has only NEA matching the reference. Each set must also meet a minimum count.

100% consistency is not sufficient: it is blind to ambiguous indels, whose alleles both
match the reference, so their orientation is untested. Trust therefore also requires that
few of the checkable ambiguous indels are suspicious -- their EAF matches the complement
of the opposite orientation's panel frequency (decide_ambiguous_indels resolves them to a
swap). A table without an EAF column cannot run that check at all and is refused trust.
"""

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.genome_reference_harmonization.allele_classes import (
    ALLELE_CLASS_COL,
    CLASS_INDEL_BOTH,
    CLASS_INDEL_EA_ONLY,
    CLASS_INDEL_NEA_ONLY,
    CLASS_INDEL_NOT_ON_REFERENCE,
    CLASS_NEA_REF,
    IS_PALINDROMIC_SNV_COL,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.ambiguous_indels import (
    INDEL_ACTION_COL,
    INDEL_DROP_REASON_COL,
    decide_ambiguous_indels,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.outcomes import (
    ACTION_SWAP,
    DROP_AMBIGUOUS_INDEL_AF_INDECISIVE,
    DROP_AMBIGUOUS_INDEL_AF_MISMATCH,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

_CONSISTENT_SNVS = "consistent_snvs"
_INCONSISTENT_SNVS = "inconsistent_snvs"
_CONSISTENT_INDELS = "consistent_indels"
_INCONSISTENT_INDELS = "inconsistent_indels"
_SUSPICIOUS = "suspicious"
_CHECKABLE = "checkable"


@frozen(slots=True)
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


@frozen(slots=True)
class SuspiciousIndelCounts:
    suspicious: (
        int  # ambiguous indels whose panel frequency contradicts the source orientation
    )
    checkable: int  # ambiguous indels with EAF and at least one panel record

    def __add__(self, other: "SuspiciousIndelCounts") -> "SuspiciousIndelCounts":
        return SuspiciousIndelCounts(
            suspicious=self.suspicious + other.suspicious,
            checkable=self.checkable + other.checkable,
        )

    @classmethod
    def zero(cls) -> "SuspiciousIndelCounts":
        return cls(suspicious=0, checkable=0)

    @property
    def fraction(self) -> float:
        return self.suspicious / self.checkable if self.checkable else 0.0


@frozen(slots=True)
class TrustEvidence:
    counts: TrustCounts
    suspicious: SuspiciousIndelCounts
    eaf_present: bool


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


def count_suspicious_indels(
    classified: pl.DataFrame,
    panel: pl.DataFrame,
    options: GenomeReferenceHarmonizationOptions,
) -> SuspiciousIndelCounts:
    """Suspicious and checkable ambiguous indels in a classify_alleles frame that carries EAF."""
    ambiguous = classified.filter(pl.col(ALLELE_CLASS_COL) == CLASS_INDEL_BOTH)
    if ambiguous.height == 0 or GWASLAB_EFFECT_ALLELE_FREQ_COL not in ambiguous.columns:
        return SuspiciousIndelCounts.zero()
    rows = ambiguous.select(
        GWASLAB_POS_COL,
        GWASLAB_EFFECT_ALLELE_COL,
        GWASLAB_NON_EFFECT_ALLELE_COL,
        pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL).cast(pl.Float64),
    )
    decisions = decide_ambiguous_indels(rows, panel, options)
    reason = pl.col(INDEL_DROP_REASON_COL)
    # Checkable = had EAF and a panel record: kept, swapped, or dropped for a frequency
    # reason (mismatch/indecisive), but not dropped for NO_EAF or NOT_IN_PANEL.
    checkable = reason.is_null() | reason.is_in(
        [DROP_AMBIGUOUS_INDEL_AF_MISMATCH, DROP_AMBIGUOUS_INDEL_AF_INDECISIVE]
    )
    summary = decisions.select(
        (pl.col(INDEL_ACTION_COL) == ACTION_SWAP).sum().alias(_SUSPICIOUS),
        checkable.sum().alias(_CHECKABLE),
    ).row(0, named=True)
    return SuspiciousIndelCounts(
        suspicious=int(summary[_SUSPICIOUS]), checkable=int(summary[_CHECKABLE])
    )


def decide_trust(
    evidence: TrustEvidence, options: GenomeReferenceHarmonizationOptions
) -> bool:
    counts = evidence.counts
    consistent = (
        counts.inconsistent_snvs == 0
        and counts.inconsistent_indels == 0
        and counts.consistent_snvs >= options.min_checkable_snvs
        and counts.consistent_indels >= options.min_checkable_indels
    )
    if not consistent:
        return False
    if not evidence.eaf_present:
        # No EAF column: ambiguous-indel orientation cannot be assessed, so refuse trust
        # (decision, 2026-09-16). Untrusted resolution then drops them (NO_EAF).
        return False
    if evidence.suspicious.checkable < options.min_checkable_ambiguous_indels:
        return True
    return evidence.suspicious.fraction <= options.max_suspicious_indel_fraction
