"""
Inspect DecodeME variants that 1000 Genomes "eur" and Pan-UKBB "ukb_eur" both keep, but in
opposite orientations.

For each conflict, prints: the input row (EA, NEA, EAF, BETA), each panel's records at the
position with the frequency the harmonizer used, the hg19 reference around the position, and
the Pan-UKBB manifest's own QC and gnomAD-genomes frequency for the site (an independent
third opinion, since the manifest carries gnomad_genomes_af_EUR).

Usage:
    pixi r python -m experiments.claude.gnomad_af_reference.inspect_orientation_conflicts \
        2>&1 | tee experiments/claude/gnomad_af_reference/inspect_orientation_conflicts.log
"""

import polars as pl
import pysam

from experiments.claude.gnomad_af_reference.decode_me_panel_comparison import (
    KEPT,
    OUTCOME_COL,
    PanelChoice,
    resolve_all,
    select_choices,
)
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.pan_ukbb.pan_ukbb_variant_manifest import (
    PAN_UKBB_VARIANT_MANIFEST,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
    gwaslab_code_to_contig_name,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    PanelTable,
    resolve_panel_af_col,
    scan_sumstats_as_polars,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.pan_ukbb.pan_ukbb_allele_frequency_panel_task import (
    scan_manifest,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SNPID_COL,
)

LABELS = ["1000g_eur", "pan_ukbb_eur"]
MANIFEST_COLUMNS = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "info",
    "high_quality",
    "af_EUR",
    "gnomad_genomes_af_EUR",
]
_FLANK = 12


def _file_path(asset: object) -> FileAsset:
    assert isinstance(asset, FileAsset)
    return asset


def panel_rows(
    choice: PanelChoice, panel: PanelTable, sites: pl.DataFrame
) -> pl.DataFrame:
    return (
        pl.scan_parquet(panel.path)
        .select(
            GWASLAB_CHROM_COL,
            pl.col(GWASLAB_POS_COL).cast(pl.Int64),
            PANEL_REF_COL,
            PANEL_ALT_COL,
            pl.col(panel.af_col).alias("panel_af"),
        )
        .join(sites.lazy(), on=[GWASLAB_CHROM_COL, GWASLAB_POS_COL], how="semi")
        .with_columns(pl.lit(choice.label).alias("panel"))
        .collect()
    )


def main() -> None:
    choices = select_choices(LABELS)
    pre_task = (
        DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.pre_harmonization_table_task
    )
    assets = DEFAULT_RUNNER.run(
        [
            pre_task,
            UCSC_HG19_INDEXED_FASTA,
            PAN_UKBB_VARIANT_MANIFEST,
            *[choice.panel_task for choice in choices],
        ]
    )
    fasta_asset = assets[UCSC_HG19_INDEXED_FASTA.asset_id]
    assert isinstance(fasta_asset, DirectoryAsset)
    fasta = IndexedFasta.open(fasta_asset.path)
    sumstats = scan_sumstats_as_polars(
        assets[pre_task.asset_id], pre_task.meta, IdentityPipe()
    )
    panels = {
        choice.label: PanelTable(
            path=_file_path(assets[choice.panel_task.asset_id]).path,
            af_col=resolve_panel_af_col(choice.panel_task, choice.ancestry),
        )
        for choice in choices
    }
    outcomes = {
        choice.label: resolve_all(choice, sumstats, fasta, panels[choice.label])
        for choice in choices
    }
    left, right = LABELS
    conflicts = (
        outcomes[left]
        .select(GWASLAB_SNPID_COL, GWASLAB_CHROM_COL, GWASLAB_POS_COL, OUTCOME_COL)
        .join(
            outcomes[right].select(GWASLAB_SNPID_COL, OUTCOME_COL),
            on=GWASLAB_SNPID_COL,
            suffix="_" + right,
        )
        .filter(
            pl.col(OUTCOME_COL).str.starts_with(KEPT)
            & pl.col(OUTCOME_COL + "_" + right).str.starts_with(KEPT)
            & (pl.col(OUTCOME_COL) != pl.col(OUTCOME_COL + "_" + right))
        )
        .rename({OUTCOME_COL: OUTCOME_COL + "_" + left})
        .sort(GWASLAB_CHROM_COL, GWASLAB_POS_COL)
    )
    print(f"conflicts: {conflicts.height}")
    sites = conflicts.select(
        pl.col(GWASLAB_CHROM_COL).cast(pl.Int32), pl.col(GWASLAB_POS_COL).cast(pl.Int64)
    ).unique()
    inputs = (
        sumstats.filter(
            pl.col(GWASLAB_SNPID_COL).is_in(conflicts[GWASLAB_SNPID_COL].implode())
        )
        .select(
            GWASLAB_SNPID_COL,
            GWASLAB_EFFECT_ALLELE_COL,
            GWASLAB_NON_EFFECT_ALLELE_COL,
            GWASLAB_EFFECT_ALLELE_FREQ_COL,
            GWASLAB_BETA_COL,
        )
        .collect()
    )
    records = pl.concat(
        [panel_rows(choice, panels[choice.label], sites) for choice in choices]
    )
    manifest = (
        scan_manifest(
            assets[PAN_UKBB_VARIANT_MANIFEST.asset_id], PAN_UKBB_VARIANT_MANIFEST.meta
        )
        .select(MANIFEST_COLUMNS)
        .join(
            sites.lazy().select(
                pl.col(GWASLAB_CHROM_COL)
                .map_elements(gwaslab_code_to_contig_name, return_dtype=pl.String)
                .alias("chrom"),
                pl.col(GWASLAB_POS_COL).alias("pos"),
            ),
            on=["chrom", "pos"],
            how="semi",
        )
        .collect(engine="streaming")
    )
    genome = pysam.FastaFile(str(fasta.fasta_path))
    with pl.Config(tbl_rows=20, tbl_cols=12, fmt_str_lengths=30, tbl_width_chars=160):
        for row in conflicts.join(inputs, on=GWASLAB_SNPID_COL).iter_rows(named=True):
            chrom, pos = row[GWASLAB_CHROM_COL], row[GWASLAB_POS_COL]
            contig = "chr" + gwaslab_code_to_contig_name(chrom)
            context = genome.fetch(contig, pos - 1 - _FLANK, pos - 1 + _FLANK)
            print(
                f"\n### {row[GWASLAB_SNPID_COL]} {chrom}:{pos} "
                f"input EA={row[GWASLAB_EFFECT_ALLELE_COL]} "
                f"NEA={row[GWASLAB_NON_EFFECT_ALLELE_COL]} "
                f"EAF={row[GWASLAB_EFFECT_ALLELE_FREQ_COL]} BETA={row[GWASLAB_BETA_COL]}"
            )
            print(
                f"{left}: {row[OUTCOME_COL + '_' + left]}   "
                f"{right}: {row[OUTCOME_COL + '_' + right]}"
            )
            print(
                f"hg19 {contig}:{pos - _FLANK}-{pos + _FLANK - 1}: "
                f"{context[:_FLANK]}[{context[_FLANK]}]{context[_FLANK + 1 :]}"
            )
            print(
                records.filter(
                    (pl.col(GWASLAB_CHROM_COL) == chrom)
                    & (pl.col(GWASLAB_POS_COL) == pos)
                )
            )
            print(
                manifest.filter(
                    (pl.col("chrom") == gwaslab_code_to_contig_name(chrom))
                    & (pl.col("pos") == pos)
                )
            )


if __name__ == "__main__":
    main()
