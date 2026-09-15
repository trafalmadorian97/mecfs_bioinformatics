"""
Inspect the ambiguous indels that V3 resolves to "swap" on build-38 DecodeME, where the
source orientation is taken as truth.

For each wrong choice at the strictest grid cell, print the row, both panel readings,
the distances, and every panel record at that position, so we can tell whether the
rule, the panel or the truth set is at fault.

Run:
  pixi r python -m experiments.claude.genome_reference_harmonization.inspect_v3_wrong_choices \
    2>&1 | tee experiments/claude/genome_reference_harmonization/inspect_v3_wrong_choices.log
"""

import attrs
import polars as pl

from experiments.claude.genome_reference_harmonization.tune_ambiguous_indel_rules import (
    BASE_OPTIONS,
    BUILD_38_TABLE,
    _ambiguous_indels,
)
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.genome_reference_harmonization.ambiguous_indels import (
    INDEL_ACTION_COL,
    decide_ambiguous_indels,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    chromosomes_to_harmonize,
    scan_sumstats_as_polars,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.outcomes import (
    ACTION_SWAP,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

STRICT_DISTANCE = 0.02
STRICT_MARGIN = 0.3
PANEL_RECORDS_COL = "panel_records_at_position"


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
    options = attrs.evolve(
        BASE_OPTIONS,
        indel_max_af_distance=STRICT_DISTANCE,
        indel_min_af_margin=STRICT_MARGIN,
    )
    wrong_parts = []
    for chrom in chromosomes_to_harmonize(sumstats, fasta, BASE_OPTIONS):
        item = _ambiguous_indels(sumstats, chrom, fasta, panel_asset.path)
        decisions = decide_ambiguous_indels(item.rows, item.panel, options)
        wrong = item.rows.filter(decisions[INDEL_ACTION_COL] == ACTION_SWAP)
        if wrong.height == 0:
            continue
        records = (
            item.panel.join(wrong.select(GWASLAB_POS_COL).unique(), on=GWASLAB_POS_COL)
            .group_by(GWASLAB_POS_COL)
            .agg(
                pl.format(
                    "{}>{}:{}",
                    pl.col(PANEL_REF_COL),
                    pl.col(PANEL_ALT_COL),
                    pl.col(PANEL_AF_COL).round(3),
                )
                .str.join("; ")
                .alias(PANEL_RECORDS_COL)
            )
        )
        wrong_parts.append(
            wrong.with_columns(pl.lit(chrom).alias(GWASLAB_CHROM_COL)).join(
                records, on=GWASLAB_POS_COL, how="left"
            )
        )
    table = pl.concat(wrong_parts).select(
        GWASLAB_CHROM_COL,
        GWASLAB_POS_COL,
        GWASLAB_EFFECT_ALLELE_COL,
        GWASLAB_NON_EFFECT_ALLELE_COL,
        GWASLAB_EFFECT_ALLELE_FREQ_COL,
        PANEL_RECORDS_COL,
    )
    with pl.Config(tbl_rows=100, tbl_cols=-1, tbl_width_chars=250, fmt_str_lengths=120):
        print(f"wrong choices at distance {STRICT_DISTANCE}, margin {STRICT_MARGIN}:")
        print(table)


if __name__ == "__main__":
    main()
