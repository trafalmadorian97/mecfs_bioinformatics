"""
Variants of DecodeME that two panel choices both keep, but in opposite orientations.

Harmonizes DecodeME under the four choices of decode_me_panel_comparison, and for every
variant kept in opposite orientations by at least one pair of choices prints: the vote
pattern (which choices keep which orientation), the input row, each panel's records at the
position, and whether the Pan-UKBB manifest row carries the known swap signature (chr21,
chr22 or X; info 1.0; no gnomAD-genomes frequency). A summary of vote patterns comes first.

Usage:
    pixi r python -m experiments.claude.gnomad_af_reference.inspect_four_way_conflicts \
        2>&1 | tee experiments/claude/gnomad_af_reference/inspect_four_way_conflicts.log
"""

import polars as pl

from experiments.claude.gnomad_af_reference.decode_me_panel_comparison import (
    KEPT,
    OUTCOME_COL,
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
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SNPID_COL,
)

SWAP_SIGNATURE_COL = "pan_ukbb_swap_signature"
VOTES_COL = "votes"
_SWAP_CHROMOSOMES = ["21", "22", "X"]


def _file(asset: object) -> FileAsset:
    assert isinstance(asset, FileAsset)
    return asset


def kept_orientation(column: str) -> pl.Expr:
    """The kept EA/NEA, or null when the choice dropped the variant."""
    return (
        pl.when(pl.col(column).str.starts_with(KEPT))
        .then(pl.col(column).str.strip_prefix(KEPT + ":"))
        .otherwise(None)
    )


def main() -> None:
    choices = select_choices([])
    labels = [choice.label for choice in choices]
    pre_task = (
        DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.pre_harmonization_table_task
    )
    panel_tasks = list(
        {choice.panel_task.asset_id: choice.panel_task for choice in choices}.values()
    )
    assets = DEFAULT_RUNNER.run(
        [pre_task, UCSC_HG19_INDEXED_FASTA, PAN_UKBB_VARIANT_MANIFEST, *panel_tasks]
    )
    fasta_asset = assets[UCSC_HG19_INDEXED_FASTA.asset_id]
    assert isinstance(fasta_asset, DirectoryAsset)
    fasta = IndexedFasta.open(fasta_asset.path)
    sumstats = scan_sumstats_as_polars(
        assets[pre_task.asset_id], pre_task.meta, IdentityPipe()
    )
    panels = {
        choice.label: PanelTable(
            path=_file(assets[choice.panel_task.asset_id]).path,
            af_col=resolve_panel_af_col(choice.panel_task, choice.ancestry),
        )
        for choice in choices
    }
    joined: pl.DataFrame | None = None
    for choice in choices:
        outcome = resolve_all(choice, sumstats, fasta, panels[choice.label]).select(
            GWASLAB_SNPID_COL,
            GWASLAB_CHROM_COL,
            GWASLAB_POS_COL,
            kept_orientation(OUTCOME_COL).alias(choice.label),
        )
        joined = (
            outcome
            if joined is None
            else joined.join(
                outcome.drop(GWASLAB_CHROM_COL, GWASLAB_POS_COL), on=GWASLAB_SNPID_COL
            )
        )
    assert joined is not None
    kept_values = pl.concat_list([pl.col(label) for label in labels]).list.drop_nulls()
    conflicts = joined.filter(kept_values.list.n_unique() > 1).with_columns(
        pl.concat_str(
            [pl.lit(label + "=") + pl.col(label).fill_null("drop") for label in labels],
            separator="  ",
        ).alias(VOTES_COL)
    )
    sites = conflicts.select(
        pl.col(GWASLAB_CHROM_COL)
        .map_elements(gwaslab_code_to_contig_name, return_dtype=pl.String)
        .alias("chrom"),
        pl.col(GWASLAB_POS_COL).alias("pos"),
    ).unique()
    signature = (
        scan_manifest(
            assets[PAN_UKBB_VARIANT_MANIFEST.asset_id], PAN_UKBB_VARIANT_MANIFEST.meta
        )
        .filter(
            pl.col("chrom").is_in(_SWAP_CHROMOSOMES)
            & (pl.col("info") == 1.0)
            & pl.col("gnomad_genomes_af_EUR").is_null()
        )
        .select("chrom", "pos")
        .join(sites.lazy(), on=["chrom", "pos"], how="semi")
        .unique()
        .with_columns(pl.lit(True).alias(SWAP_SIGNATURE_COL))
        .collect(engine="streaming")
    )
    conflicts = (
        conflicts.with_columns(
            pl.col(GWASLAB_CHROM_COL)
            .map_elements(gwaslab_code_to_contig_name, return_dtype=pl.String)
            .alias("chrom"),
            pl.col(GWASLAB_POS_COL).alias("pos"),
        )
        .join(signature, on=["chrom", "pos"], how="left")
        .with_columns(pl.col(SWAP_SIGNATURE_COL).fill_null(False))
        .drop("chrom", "pos")
        .sort(GWASLAB_CHROM_COL, GWASLAB_POS_COL)
    )
    print(
        f"variants kept in opposite orientations by some pair of choices: {conflicts.height}"
    )
    with pl.Config(tbl_rows=60, fmt_str_lengths=200, tbl_width_chars=250):
        print(
            "\n== who sides with whom (orientation agreement pattern) x Pan-UKBB swap signature"
        )
        pattern = pl.concat_str(
            [
                pl.when(pl.col(label).is_null())
                .then(pl.lit("-"))
                .when(pl.col(label) == pl.coalesce([pl.col(name) for name in labels]))
                .then(pl.lit("A"))
                .otherwise(pl.lit("B"))
                for label in labels
            ]
        ).alias("pattern_" + "/".join(labels))
        print(
            conflicts.group_by(pattern, SWAP_SIGNATURE_COL)
            .len()
            .sort("len", descending=True)
        )
        inputs = (
            sumstats.filter(
                pl.col(GWASLAB_SNPID_COL).is_in(conflicts[GWASLAB_SNPID_COL].implode())
            )
            .select(
                GWASLAB_SNPID_COL,
                GWASLAB_EFFECT_ALLELE_COL,
                GWASLAB_NON_EFFECT_ALLELE_COL,
                GWASLAB_EFFECT_ALLELE_FREQ_COL,
            )
            .collect()
        )
        conflict_positions = conflicts[GWASLAB_POS_COL].unique().to_list()
        records = pl.concat(
            [
                pl.scan_parquet(panels[label].path)
                # A pushed-down position filter keeps the 224M-row gnomAD scan small.
                .filter(pl.col(GWASLAB_POS_COL).is_in(conflict_positions))
                .select(
                    GWASLAB_CHROM_COL,
                    pl.col(GWASLAB_POS_COL).cast(pl.Int64),
                    pl.lit(label).alias("panel"),
                    PANEL_REF_COL,
                    PANEL_ALT_COL,
                    pl.col(panels[label].af_col).alias("panel_af"),
                )
                .join(
                    conflicts.lazy().select(
                        pl.col(GWASLAB_CHROM_COL).cast(pl.Int32),
                        pl.col(GWASLAB_POS_COL).cast(pl.Int64),
                    ),
                    on=[GWASLAB_CHROM_COL, GWASLAB_POS_COL],
                    how="semi",
                )
                .collect()
                for label in labels
            ]
        )
        multiallelic = (
            records.filter(pl.col("panel") == "gnomad_nfe")
            .group_by(GWASLAB_CHROM_COL, GWASLAB_POS_COL)
            .len("gnomad_records_at_position")
        )
        print(
            "\n== pattern x swap signature x gnomAD records at the position (>1: multiallelic)"
        )
        print(
            conflicts.with_columns(pl.col(GWASLAB_CHROM_COL).cast(pl.Int32))
            .join(multiallelic, on=[GWASLAB_CHROM_COL, GWASLAB_POS_COL], how="left")
            .group_by(
                pattern,
                SWAP_SIGNATURE_COL,
                (pl.col("gnomad_records_at_position").fill_null(0) > 1).alias(
                    "gnomad_multiallelic"
                ),
            )
            .len()
            .sort("len", descending=True)
        )
        print("\n== conflicts without the Pan-UKBB swap signature, in detail")
        for row in (
            conflicts.filter(~pl.col(SWAP_SIGNATURE_COL))
            .join(inputs, on=GWASLAB_SNPID_COL)
            .iter_rows(named=True)
        ):
            chrom, pos = row[GWASLAB_CHROM_COL], row[GWASLAB_POS_COL]
            print(
                f"\n### {row[GWASLAB_SNPID_COL]} input EA={row[GWASLAB_EFFECT_ALLELE_COL]} "
                f"NEA={row[GWASLAB_NON_EFFECT_ALLELE_COL]} "
                f"EAF={row[GWASLAB_EFFECT_ALLELE_FREQ_COL]:.4f}"
            )
            print(row[VOTES_COL])
            print(
                records.filter(
                    (pl.col(GWASLAB_CHROM_COL) == chrom)
                    & (pl.col(GWASLAB_POS_COL) == pos)
                )
            )


if __name__ == "__main__":
    main()
