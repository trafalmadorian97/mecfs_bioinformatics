"""
V3: tune indel_max_af_distance and indel_min_af_margin on a truth set.

DecodeME build 38 is 100% reference-consistent against hg38, so each ambiguous indel's
source orientation is the truth. Every ambiguous indel is resolved with the stringent
rules as if the table were untrusted, over a grid of options. Choosing "swap" is a wrong
choice. The chosen defaults are the grid cell with zero wrong choices that keeps the most
indels; ties go to the larger margin, then the smaller distance.

Run:
  pixi r python -m experiments.claude.genome_reference_harmonization.tune_ambiguous_indel_rules \
    2>&1 | tee experiments/claude/genome_reference_harmonization/tune_ambiguous_indel_rules.log
"""

from pathlib import Path

import attrs
import polars as pl
from attrs import frozen

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_sumstats_minimal_processing import (
    DECODE_ME_GWAS_1_SUMSTATS_MINIMAL_FILTERING,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.genome_reference_harmonization.allele_classes import (
    ALLELE_CLASS_COL,
    CLASS_INDEL_BOTH,
    classify_alleles,
    prepare_alleles,
    valid_alleles_expr,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.ambiguous_indels import (
    INDEL_ACTION_COL,
    INDEL_DROP_REASON_COL,
    decide_ambiguous_indels,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    ParquetPanelLoader,
    chromosomes_to_harmonize,
    count_trust_evidence_genome_wide,
    scan_sumstats_as_polars,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.outcomes import (
    ACTION_KEEP,
    ACTION_SWAP,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.trust import (
    decide_trust,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

GRID_DISTANCES = (0.02, 0.05, 0.1, 0.15, 0.2)
GRID_MARGINS = (0.02, 0.05, 0.1, 0.2, 0.3)
BASE_OPTIONS = GenomeReferenceHarmonizationOptions()
ROW_COLUMNS = [
    GWASLAB_POS_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
]

BUILD_38_TABLE = GwasLabSumstatsToTableTask.create_from_source_task(
    source_tsk=DECODE_ME_GWAS_1_SUMSTATS_MINIMAL_FILTERING,
    asset_id="experiment_decode_me_gwas_1_build_38_table",
    sub_dir="processed",
)


@frozen
class ChromosomeAmbiguousIndels:
    rows: pl.DataFrame
    panel: pl.DataFrame


def _ambiguous_indels(
    sumstats: pl.LazyFrame, chrom: int, fasta: IndexedFasta, panel_path: Path
) -> ChromosomeAmbiguousIndels:
    rows = (
        sumstats.filter(pl.col(GWASLAB_CHROM_COL) == chrom)
        .select(ROW_COLUMNS)
        .collect(engine="streaming")
    )
    classified = classify_alleles(
        prepare_alleles(rows).filter(valid_alleles_expr()),
        fasta=fasta,
        chrom=chrom,
        max_gather_bytes=BASE_OPTIONS.max_gather_bytes,
    )
    ambiguous = classified.filter(pl.col(ALLELE_CLASS_COL) == CLASS_INDEL_BOTH).select(
        ROW_COLUMNS
    )
    panel = ParquetPanelLoader(panel_path=panel_path, chrom=chrom)(
        ambiguous[GWASLAB_POS_COL]
    )
    return ChromosomeAmbiguousIndels(rows=ambiguous, panel=panel)


def main() -> None:
    assets = DEFAULT_RUNNER.run(
        [
            BUILD_38_TABLE,
            UCSC_HG38_INDEXED_FASTA,
            THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
        ]
    )
    fasta_asset = assets[UCSC_HG38_INDEXED_FASTA.asset_id]
    panel_asset = assets[THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES.asset_id]
    assert isinstance(fasta_asset, DirectoryAsset) and isinstance(
        panel_asset, FileAsset
    )
    fasta = IndexedFasta.open(fasta_asset.path)
    sumstats = scan_sumstats_as_polars(
        assets[BUILD_38_TABLE.asset_id], BUILD_38_TABLE.meta, IdentityPipe()
    )
    chromosomes = chromosomes_to_harmonize(sumstats, fasta, BASE_OPTIONS)
    evidence = count_trust_evidence_genome_wide(
        sumstats, chromosomes, fasta, panel_asset.path, BASE_OPTIONS
    )
    print(f"build-38 trust counts: {evidence.counts}")
    assert decide_trust(evidence, BASE_OPTIONS), (
        "build 38 is not trusted, so it is not a truth set"
    )
    per_chromosome = [
        _ambiguous_indels(sumstats, chrom, fasta, panel_asset.path)
        for chrom in chromosomes
    ]
    n_ambiguous = sum(item.rows.height for item in per_chromosome)
    print(f"ambiguous indels: {n_ambiguous:,}")
    results = []
    for distance in GRID_DISTANCES:
        for margin in GRID_MARGINS:
            options = attrs.evolve(
                BASE_OPTIONS, indel_max_af_distance=distance, indel_min_af_margin=margin
            )
            decisions = pl.concat(
                [
                    decide_ambiguous_indels(item.rows, item.panel, options)
                    for item in per_chromosome
                ]
            )
            reasons = dict(decisions.group_by(INDEL_DROP_REASON_COL).len().rows())
            results.append(
                {
                    "distance": distance,
                    "margin": margin,
                    "kept": int((decisions[INDEL_ACTION_COL] == ACTION_KEEP).sum()),
                    "wrong": int((decisions[INDEL_ACTION_COL] == ACTION_SWAP).sum()),
                    **{
                        f"drop_{reason}": count
                        for reason, count in reasons.items()
                        if reason is not None
                    },
                }
            )
    table = pl.DataFrame(results).fill_null(0)
    with pl.Config(tbl_rows=40, tbl_cols=-1, tbl_width_chars=250):
        print(table)
    candidates = table.filter(pl.col("wrong") == 0).sort(
        ["kept", "margin", "distance"], descending=[True, True, False]
    )
    assert candidates.height > 0, (
        "no grid cell has zero wrong choices; report to the user"
    )
    print("chosen defaults:", candidates.row(0, named=True))


if __name__ == "__main__":
    main()
