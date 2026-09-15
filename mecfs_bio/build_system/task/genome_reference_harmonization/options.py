"""Options for genome-reference harmonization, validated at construction."""

from attrs import frozen

from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    DEFAULT_MAX_GATHER_BYTES,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.flip import (
    ExtraColumnRule,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_CODE_FOR_NAME


@frozen
class GenomeReferenceHarmonizationOptions:
    """
    min_checkable_snvs, min_checkable_indels: a table is trusted only if it has at least
        this many checkable SNVs (non-palindromic) and indels, all consistent with NEA as the
        reference allele.
    palindrome_maf_threshold, panel_maf_threshold: an untrusted palindromic SNV is resolved
        by frequency only when both its summary-statistic MAF and its panel MAF are at most
        these (gwaslab's defaults).
    indel_max_af_distance, indel_min_af_margin: stringent rules for untrusted ambiguous indels.
        A reading is chosen only if its predicted EAF is within the distance and beats the
        other reading by at least the margin.
    keep_unresolved_palindromes: keep untrusted palindromic SNVs whose strand cannot be
        resolved, instead of dropping them. Never applies to indels.
    excluded_chromosomes: gwaslab chromosome codes whose rows are dropped (MT by default,
        since UCSC hg19 chrM is not rCRS).
    extra_column_rules: flip rules for dataset-specific columns.
    max_gather_bytes: memory budget for one reference-gather slice.
    """

    min_checkable_snvs: int = 10_000
    min_checkable_indels: int = 1_000
    palindrome_maf_threshold: float = 0.4
    panel_maf_threshold: float = 0.4
    indel_max_af_distance: float = 0.1
    indel_min_af_margin: float = 0.2
    keep_unresolved_palindromes: bool = False
    excluded_chromosomes: tuple[int, ...] = (GWASLAB_CHROM_CODE_FOR_NAME["MT"],)
    extra_column_rules: tuple[ExtraColumnRule, ...] = ()
    max_gather_bytes: int = DEFAULT_MAX_GATHER_BYTES

    def __attrs_post_init__(self) -> None:
        assert self.min_checkable_snvs >= 1 and self.min_checkable_indels >= 1
        assert 0 < self.palindrome_maf_threshold < 0.5
        assert 0 < self.panel_maf_threshold < 0.5
        assert 0 < self.indel_max_af_distance < 1
        assert 0 < self.indel_min_af_margin < 1, (
            "a strictly positive margin keeps the keep and flip readings from both being chosen"
        )
        assert self.max_gather_bytes > 0
