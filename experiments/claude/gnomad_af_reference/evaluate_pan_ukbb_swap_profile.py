"""
Per chromosome: how often does a Pan-UKBB manifest row's REF/ALT orientation disagree with
gnomAD v2.1.1 genomes and with 1000 Genomes EUR?

The proposed Pan-UKBB fix drops rows with a "swap profile": genotyped (info == 1.0) and no
gnomAD-genomes frequency in the manifest, on chr21, chr22 and X only. Restricting a rule to
three chromosomes needs evidence, so this script measures the disagreement rate on every
chromosome, for the profile rows and for genotyped rows that do have a gnomAD frequency
(the control).

For each manifest row and each reference panel, records at the same position are looked up:
- agree: the panel has REF=ref, ALT=alt (and not the reverse);
- disagree: the panel has REF=alt, ALT=ref (and not the forward record);
- both: the panel has both records (two distinct variants, e.g. an insertion and a deletion);
- absent: neither.
frac_disagree is disagree / (agree + disagree). For disagreeing rows, the frequency check
says whether af_EUR tracks 1 - panel AF (a label swap) rather than panel AF.

Usage:
    pixi r python -m experiments.claude.gnomad_af_reference.evaluate_pan_ukbb_swap_profile \
        2>&1 | tee experiments/claude/gnomad_af_reference/evaluate_pan_ukbb_swap_profile.log
"""

from pathlib import Path

import polars as pl
from attrs import frozen

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.gnomad.gnomad_allele_frequency_panels import (
    GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES,
)
from mecfs_bio.assets.reference_data.pan_ukbb.pan_ukbb_variant_manifest import (
    PAN_UKBB_VARIANT_MANIFEST,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.genome_reference_harmonization.pan_ukbb.pan_ukbb_allele_frequency_panel_task import (
    scan_manifest,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL

CHROMOSOMES = [str(number) for number in range(1, 23)] + ["X"]
GROUP_COL = "group"
PROFILE = "profile: info 1, no gnomAD AF"
CONTROL = "control: info 1, gnomAD AF"
AGREE = "agree"
DISAGREE = "disagree"
BOTH = "both"
ABSENT = "absent"


@frozen(slots=True)
class ReferencePanel:
    label: str
    path: Path
    af_col: str


def _file(asset: Asset) -> FileAsset:
    assert isinstance(asset, FileAsset)
    return asset


def chrom_code(chrom: str) -> int:
    return 23 if chrom == "X" else int(chrom)


def genotyped_rows(manifest: pl.LazyFrame) -> pl.DataFrame:
    """Genotyped manifest rows, labelled profile or control."""
    return (
        manifest.filter(pl.col("info") == 1.0)
        .select(
            "chrom",
            "pos",
            "ref",
            "alt",
            "af_EUR",
            pl.when(pl.col("gnomad_genomes_af_EUR").is_null())
            .then(pl.lit(PROFILE))
            .otherwise(pl.lit(CONTROL))
            .alias(GROUP_COL),
            pl.when(
                (pl.col("ref").str.len_bytes() == 1)
                & (pl.col("alt").str.len_bytes() == 1)
            )
            .then(pl.lit("snv"))
            .otherwise(pl.lit("indel"))
            .alias("allele_class"),
        )
        .collect(engine="streaming")
    )


def orientation(rows: pl.DataFrame, panel: ReferencePanel, chrom: str) -> pl.DataFrame:
    """Per row: the panel's orientation verdict and its AF in the manifest's orientation."""
    records = (
        pl.scan_parquet(panel.path)
        .filter(pl.col(GWASLAB_CHROM_COL) == chrom_code(chrom))
        .filter(pl.col(GWASLAB_POS_COL).is_in(rows["pos"].unique().implode()))
        .select(
            pl.col(GWASLAB_POS_COL).cast(pl.Int64).alias("pos"),
            PANEL_REF_COL,
            PANEL_ALT_COL,
            pl.col(panel.af_col).cast(pl.Float64).alias("panel_af"),
        )
        .collect()
    )
    forward = rows.join(
        records.rename({PANEL_REF_COL: "ref", PANEL_ALT_COL: "alt"}),
        on=["pos", "ref", "alt"],
    ).select("row", pl.col("panel_af").alias("forward_af"))
    reverse = rows.join(
        records.rename({PANEL_REF_COL: "alt", PANEL_ALT_COL: "ref"}),
        on=["pos", "ref", "alt"],
    ).select("row", pl.col("panel_af").alias("reverse_af"))
    joined = rows.join(forward, on="row", how="left").join(
        reverse, on="row", how="left"
    )
    assert joined.height == rows.height, f"{panel.label}: duplicate panel records"
    has_forward = pl.col("forward_af").is_not_null()
    has_reverse = pl.col("reverse_af").is_not_null()
    return joined.select(
        "row",
        pl.when(has_forward & has_reverse)
        .then(pl.lit(BOTH))
        .when(has_forward)
        .then(pl.lit(AGREE))
        .when(has_reverse)
        .then(pl.lit(DISAGREE))
        .otherwise(pl.lit(ABSENT))
        .alias(panel.label),
        # For a disagreeing row: does af_EUR match 1 - panel AF (a label swap)?
        (
            has_reverse
            & ~has_forward
            & (
                (pl.col("af_EUR") - (1 - pl.col("reverse_af"))).abs()
                < (pl.col("af_EUR") - pl.col("reverse_af")).abs()
            )
        ).alias(panel.label + "_af_says_swap"),
    )


def summarize(verdicts: pl.DataFrame, labels: list[str]) -> pl.DataFrame:
    aggregations = []
    for label in labels:
        agree = (pl.col(label) == AGREE).sum()
        disagree = (pl.col(label) == DISAGREE).sum()
        aggregations += [
            agree.alias(f"{label}_agree"),
            disagree.alias(f"{label}_disagree"),
            (pl.col(label) == BOTH).sum().alias(f"{label}_both"),
            (disagree / (agree + disagree)).round(3).alias(f"{label}_frac_disagree"),
            pl.col(label + "_af_says_swap").sum().alias(f"{label}_disagree_af_swap"),
        ]
    return verdicts.group_by(
        "chrom", GROUP_COL, "allele_class", maintain_order=True
    ).agg(pl.len().alias("n"), *aggregations)


def main() -> None:
    assets = DEFAULT_RUNNER.run(
        [
            PAN_UKBB_VARIANT_MANIFEST,
            GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES,
            THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
        ]
    )
    panels = [
        ReferencePanel(
            label="gnomad",
            path=_file(
                assets[GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES.asset_id]
            ).path,
            af_col="AF_nfe",
        ),
        ReferencePanel(
            label="1000g",
            path=_file(
                assets[THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES.asset_id]
            ).path,
            af_col="AF",
        ),
    ]
    labels = [panel.label for panel in panels]
    manifest = scan_manifest(
        assets[PAN_UKBB_VARIANT_MANIFEST.asset_id], PAN_UKBB_VARIANT_MANIFEST.meta
    )
    genotyped = genotyped_rows(manifest)
    parts = []
    for chrom in CHROMOSOMES:
        rows = genotyped.filter(pl.col("chrom") == chrom).with_row_index("row")
        verdicts = rows.select("row", GROUP_COL, "allele_class")
        for panel in panels:
            verdicts = verdicts.join(orientation(rows, panel, chrom), on="row")
        parts.append(
            summarize(verdicts.with_columns(pl.lit(chrom).alias("chrom")), labels)
        )
        print(f"done chr{chrom}: {rows.height} genotyped rows", flush=True)
    summary = pl.concat(parts)
    columns = ["chrom", "allele_class", "n"] + [
        f"{label}_{suffix}"
        for label in labels
        for suffix in ["agree", "disagree", "both", "frac_disagree", "disagree_af_swap"]
    ]
    with pl.Config(tbl_rows=200, tbl_cols=20, tbl_width_chars=250):
        for group in [PROFILE, CONTROL]:
            for allele_class in ["indel", "snv"]:
                print(f"\n== {group}; {allele_class}")
                print(
                    summary.filter(
                        (pl.col(GROUP_COL) == group)
                        & (pl.col("allele_class") == allele_class)
                    ).select(columns)
                )


if __name__ == "__main__":
    main()
