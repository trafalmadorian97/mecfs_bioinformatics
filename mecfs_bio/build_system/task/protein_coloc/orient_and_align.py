"""
Orient a dense protein window to the reference genome and align it exactly to the trait.

Proteomics sources are not trusted to report the reference allele as NEA (UKB-PPP's ALLELE0 is
the hg19 reference at a small fraction of sites; a CSF pQTL source reports the in-sample minor
allele). So every protein window goes through the untrusted path of the genome-reference
harmonizer: SNVs are oriented by the FASTA, while palindromic SNVs and ambiguous indels (both
orientations on the reference) are resolved by allele frequency against a panel, or dropped.

After orientation both sides have NEA = GRCh38 reference, so trait and protein are matched by an
exact join on the ordered key (CHR, POS, EA, NEA). This keeps a mirrored insertion/deletion pair
(T>TCA versus TCA>T at one position) apart, which an unordered allele key would merge.

A joined variant whose two effect allele frequencies differ by more than max_eaf_difference is
dropped: the two sides then probably do not describe the same variant (or one mislabels its
alleles). Variants with a null EAF on either side cannot be checked and are kept. The default,
DEFAULT_MAX_EAF_DIFFERENCE = 0.2, comes from 8.28M variants of a deCODE SomaScan file (Iceland)
joined to the DecodeME GWAS (UK): the population difference |EAF_a - EAF_b| has its 99.99th
percentile at 0.14 for SNVs with MAF >= 20% and falls off near 0.15; 224 variants (2.7 per
100,000, relatively more indels) exceed 0.2. Sources from the same population differ less.

With drop_multi_record_trait_positions, every trait row at a position holding more than one
trait record is dropped first. A source whose rows are all two-allele tests (deCODE keeps only
positions where Iceland has a single ALT) needs it: where the trait has several records at a
position, its row for one ALT may pool the other ALTs' carriers with REF (as bcftools norm -m-
codes them), so it does not measure the protein row's contrast. Other sources keep these rows.
"""

from collections.abc import Mapping
from pathlib import Path

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.flip import (
    resolve_column_rules,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    ParquetPanelLoader,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.resolve_chromosome import (
    DROP_REASON_COL,
    ChromosomeContext,
    resolve_chromosome,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    BETA_PROTEIN_COL,
    BETA_TRAIT_COL,
    EAF_PROTEIN_COL,
    EAF_TRAIT_COL,
    N_PROTEIN_COL,
    N_TRAIT_COL,
    SE_PROTEIN_COL,
    SE_TRAIT_COL,
    WINDOW_SCHEMA,
)

DEFAULT_MAX_EAF_DIFFERENCE = 0.2

VARIANT_KEY = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
]


@frozen(slots=True)
class OrientedWindow:
    """
    frame: the kept rows in WINDOW_SCHEMA, with NEA the reference allele.
    n_input: rows before orientation.
    drop_counts: number of dropped rows per drop reason.
    """

    frame: pl.DataFrame
    n_input: int
    drop_counts: Mapping[str, int]


def orient_window(
    rows: pl.DataFrame,
    chrom: int,
    fasta: IndexedFasta,
    panel_path: Path,
    panel_af_col: str,
    options: GenomeReferenceHarmonizationOptions,
) -> OrientedWindow:
    """Orient one chromosome's window rows with the untrusted genome-reference path."""
    assert rows.schema == WINDOW_SCHEMA, f"unexpected window schema: {rows.schema}"
    resolved = resolve_chromosome(
        rows,
        ChromosomeContext(
            chrom=chrom,
            fasta=fasta,
            trusted=False,
            rules=resolve_column_rules(
                columns=rows.columns, extra=options.extra_column_rules
            ),
            options=options,
        ),
        ParquetPanelLoader(panel_path=panel_path, af_col=panel_af_col, chrom=chrom),
    )
    dropped = resolved.filter(pl.col(DROP_REASON_COL).is_not_null())
    drop_counts = {
        reason: count
        for reason, count in dropped.group_by(DROP_REASON_COL).len().rows()
    }
    kept = (
        resolved.filter(pl.col(DROP_REASON_COL).is_null())
        .select(WINDOW_SCHEMA.names())
        .cast(WINDOW_SCHEMA)
    )
    return OrientedWindow(frame=kept, n_input=rows.height, drop_counts=drop_counts)


@frozen(slots=True)
class AlignedWindow:
    """
    frame: one row per variant present on both sides: the key CHR, POS, EA, NEA, then
        BETA/SE/EAF/N of the trait and of the protein (suffixes _TRAIT, _PROTEIN).
    n_trait, n_protein: input rows of each side.
    n_duplicate_keys_dropped: distinct keys that occurred more than once within a side,
        summed over both sides; all their rows are dropped.
    n_eaf_discordant_dropped: joined variants dropped because their EAFs differ by more than
        max_eaf_difference.
    n_trait_multi_record_dropped: trait rows dropped because their position holds more than
        one trait record (only with drop_multi_record_trait_positions).
    """

    frame: pl.DataFrame
    n_trait: int
    n_protein: int
    n_duplicate_keys_dropped: int
    n_eaf_discordant_dropped: int
    n_trait_multi_record_dropped: int


def _key_is_duplicated() -> pl.Expr:
    return pl.struct(VARIANT_KEY).is_duplicated()


def _count_duplicated_keys(frame: pl.DataFrame) -> int:
    return frame.filter(_key_is_duplicated()).select(VARIANT_KEY).unique().height


def _suffixed(frame: pl.DataFrame, names: Mapping[str, str]) -> pl.DataFrame:
    return frame.select(*VARIANT_KEY, *[pl.col(k).alias(v) for k, v in names.items()])


def align_windows(
    trait: pl.DataFrame,
    protein: pl.DataFrame,
    max_eaf_difference: float | None = DEFAULT_MAX_EAF_DIFFERENCE,
    drop_multi_record_trait_positions: bool = False,
) -> AlignedWindow:
    """
    Exact inner join of two reference-oriented windows on (CHR, POS, EA, NEA), without the
    variants whose EAFs differ by more than max_eaf_difference (None: no check), and, with
    drop_multi_record_trait_positions, without trait positions holding several records.
    """
    assert max_eaf_difference is None or 0 < max_eaf_difference < 1
    single_record_trait = (
        trait.filter(pl.len().over(GWASLAB_CHROM_COL, GWASLAB_POS_COL) == 1)
        if drop_multi_record_trait_positions
        else trait
    )
    joined = _suffixed(
        single_record_trait.filter(~_key_is_duplicated()),
        {
            GWASLAB_BETA_COL: BETA_TRAIT_COL,
            GWASLAB_SE_COL: SE_TRAIT_COL,
            GWASLAB_EFFECT_ALLELE_FREQ_COL: EAF_TRAIT_COL,
            GWASLAB_SAMPLE_SIZE_COLUMN: N_TRAIT_COL,
        },
    ).join(
        _suffixed(
            protein.filter(~_key_is_duplicated()),
            {
                GWASLAB_BETA_COL: BETA_PROTEIN_COL,
                GWASLAB_SE_COL: SE_PROTEIN_COL,
                GWASLAB_EFFECT_ALLELE_FREQ_COL: EAF_PROTEIN_COL,
                GWASLAB_SAMPLE_SIZE_COLUMN: N_PROTEIN_COL,
            },
        ),
        on=VARIANT_KEY,
        how="inner",
    )
    discordant = (
        pl.lit(False)
        if max_eaf_difference is None
        else (
            (pl.col(EAF_TRAIT_COL) - pl.col(EAF_PROTEIN_COL)).abs() > max_eaf_difference
        ).fill_null(False)
    )
    concordant = joined.filter(~discordant)
    return AlignedWindow(
        frame=concordant.sort(GWASLAB_POS_COL),
        n_trait=trait.height,
        n_protein=protein.height,
        n_duplicate_keys_dropped=_count_duplicated_keys(single_record_trait)
        + _count_duplicated_keys(protein),
        n_eaf_discordant_dropped=joined.height - concordant.height,
        n_trait_multi_record_dropped=trait.height - single_record_trait.height,
    )
