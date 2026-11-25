"""
How many rows of a deCODE SomaScan file sit at multi-row positions (a position with more than
one row: multi-allelic sites, SNVs under spanning deletions, and the "*" and X:X rows)? Reads the
local reference copy (DECODE_SMP_AP4M1_RAW_FILE).

  pixi r python experiments/claude/decode_proteomics_survey/count_multiallelic_rows.py
"""

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_index import (
    HAPMAP_3_DECODE_SMP_INDEX,
)
from mecfs_bio.assets.gwas.decode_somascan.raw.decode_reference_protein_file import (
    DECODE_SMP_AP4M1_RAW_FILE,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset

pl.Config.set_tbl_rows(20)
pl.Config.set_tbl_width_chars(200)


def main() -> None:
    assets = DEFAULT_RUNNER.run([DECODE_SMP_AP4M1_RAW_FILE, HAPMAP_3_DECODE_SMP_INDEX])
    raw = assets[DECODE_SMP_AP4M1_RAW_FILE.asset_id]
    index = assets[HAPMAP_3_DECODE_SMP_INDEX.asset_id]
    assert isinstance(raw, FileAsset) and isinstance(index, FileAsset)
    index_names = pl.read_parquet(index.path, columns=["NAME"])["NAME"]

    rows = pl.read_csv(
        raw.path,
        separator="\t",
        columns=["Chrom", "Pos", "Name", "effectAllele", "otherAllele", "ImpMAF"],
        null_values=["."],
        quote_char=None,
        schema_overrides={"Pos": pl.Int64, "ImpMAF": pl.Float64},
    )
    ea, oa = pl.col("effectAllele"), pl.col("otherAllele")
    rows = rows.with_columns(
        pl.when((ea == "*") | (oa == "*"))
        .then(pl.lit("star"))
        .when(ea == oa)
        .then(pl.lit("same_allele"))
        .when((ea.str.len_bytes() == 1) & (oa.str.len_bytes() == 1))
        .then(pl.lit("snv"))
        .otherwise(pl.lit("indel"))
        .alias("kind"),
        (pl.len().over("Chrom", "Pos") > 1).alias("multi"),
        (pl.col("ImpMAF") >= 0.01).fill_null(False).alias("common"),
        pl.col("Name").is_in(index_names.implode()).alias("in_index"),
    )

    positions = rows.select(pl.struct("Chrom", "Pos").n_unique()).item()
    multi_positions = (
        rows.filter("multi").select(pl.struct("Chrom", "Pos").n_unique()).item()
    )
    print(
        f"{rows.height} rows at {positions} positions; {multi_positions} multi-row positions"
    )
    dropped = int(rows["multi"].sum())
    print(f"rows at multi-row positions: {dropped} ({dropped / rows.height:.1%})")

    def share(mask: pl.Expr) -> pl.Expr:
        return mask.sum() / pl.len()

    by_kind = (
        rows.group_by("kind")
        .agg(
            pl.len().alias("rows"),
            pl.col("multi").sum().alias("dropped"),
            share(pl.col("multi")).alias("dropped_share"),
            pl.col("common").sum().alias("rows_maf_ge_1pct"),
            (pl.col("multi") & pl.col("common")).sum().alias("dropped_maf_ge_1pct"),
        )
        .with_columns(
            (pl.col("dropped_maf_ge_1pct") / pl.col("rows_maf_ge_1pct")).alias(
                "dropped_share_maf_ge_1pct"
            )
        )
        .sort("rows", descending=True)
    )
    print(by_kind)
    common_snvs = rows.filter(pl.col("common") & (pl.col("kind") == "snv"))
    print(
        f"common SNVs (MAF >= 1%) at multi-row positions: {int(common_snvs['multi'].sum())} "
        f"of {common_snvs.height}"
    )
    print(
        f"HapMap3 index rows at multi-row positions: "
        f"{int(rows.filter('in_index')['multi'].sum())} of {index_names.len()}"
    )


if __name__ == "__main__":
    main()
