"""
How many Pan-UKBB manifest rows are REF/ALT-swapped but pass the FASTA check?

The 295 known swaps (chr21 47, chr22 86, X 162) are genotyped sites (info 1.0) whose listed
ref is not the hg19 base. An ambiguous indel (both orientations on the FASTA) whose labels
are swapped would pass that check unnoticed. This probe profiles genotyped (info == 1.0)
rows by chromosome and allele class, and joins them to the built panel to see which survive
into it.

Usage:
    pixi r python -m experiments.claude.gnomad_af_reference.probe_undetected_pan_ukbb_swaps \
        2>&1 | tee experiments/claude/gnomad_af_reference/probe_undetected_pan_ukbb_swaps.log
"""

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.pan_ukbb.pan_ukbb_allele_frequencies import (
    PAN_UKBB_HG19_ALLELE_FREQUENCIES,
)
from mecfs_bio.assets.reference_data.pan_ukbb.pan_ukbb_variant_manifest import (
    PAN_UKBB_VARIANT_MANIFEST,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.genome_reference_harmonization.pan_ukbb.pan_ukbb_allele_frequency_panel_task import (
    scan_manifest,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL

ALLELE_CLASS = (
    pl.when((pl.col("ref").str.len_bytes() == 1) & (pl.col("alt").str.len_bytes() == 1))
    .then(pl.lit("snv"))
    .otherwise(pl.lit("indel"))
    .alias("allele_class")
)


def main() -> None:
    assets = DEFAULT_RUNNER.run(
        [PAN_UKBB_VARIANT_MANIFEST, PAN_UKBB_HG19_ALLELE_FREQUENCIES]
    )
    panel_asset = assets[PAN_UKBB_HG19_ALLELE_FREQUENCIES.asset_id]
    assert isinstance(panel_asset, FileAsset)
    genotyped = (
        scan_manifest(
            assets[PAN_UKBB_VARIANT_MANIFEST.asset_id], PAN_UKBB_VARIANT_MANIFEST.meta
        )
        .filter(pl.col("info") == 1.0)
        .select(
            "chrom",
            "pos",
            "ref",
            "alt",
            "high_quality",
            "af_EUR",
            "gnomad_genomes_af_EUR",
            ALLELE_CLASS,
        )
        .collect(engine="streaming")
    )
    panel_keys = (
        pl.scan_parquet(panel_asset.path)
        .select(
            pl.col(GWASLAB_CHROM_COL)
            .replace_strict(
                {23: "X"}, default=pl.col(GWASLAB_CHROM_COL).cast(pl.String)
            )
            .alias("chrom"),
            pl.col(GWASLAB_POS_COL).cast(pl.Int64).alias("pos"),
            pl.col("REF").alias("ref"),
            pl.col("ALT").alias("alt"),
            pl.lit(True).alias("in_panel"),
        )
        .collect(engine="streaming")
    )
    genotyped = genotyped.join(
        panel_keys, on=["chrom", "pos", "ref", "alt"], how="left"
    ).with_columns(pl.col("in_panel").fill_null(False))
    with pl.Config(tbl_rows=200):
        print(
            "== genotyped (info == 1.0) manifest rows by chromosome, class, panel membership"
        )
        print(
            genotyped.group_by("chrom", "allele_class", "in_panel", "high_quality")
            .len()
            .sort("chrom", "allele_class", "in_panel", "high_quality")
        )
        print("== genotyped rows with no gnomAD-genomes frequency, by chromosome")
        print(
            genotyped.group_by(
                "chrom", pl.col("gnomad_genomes_af_EUR").is_null().alias("no_gnomad")
            )
            .len()
            .sort("chrom", "no_gnomad")
        )


if __name__ == "__main__":
    main()
