"""
For every deCODE SNV row (reference file AP4M1): q = total frequency of the OTHER non-reference
alleles at its position, from variant_info_iceland (which lists every non-reference allele
genotyped at a position, "*" spanning deletions included). q is how much of the "rest" that the
SNV's allele is tested against is not the reference allele. Summarised for SNVs alone at their
position and for SNVs at multi-row positions, by MAF.

  pixi r python experiments/claude/decode_proteomics_survey/other_allele_mass.py
"""

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.decode_somascan.raw.decode_reference_protein_file import (
    DECODE_SMP_AP4M1_RAW_FILE,
)
from mecfs_bio.assets.gwas.decode_somascan.raw.decode_variant_info import (
    DECODE_VARIANT_INFO,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.decode_database.decode_sumstats_format import (
    chrom_code_expr,
)

pl.Config.set_tbl_rows(40)
pl.Config.set_tbl_width_chars(220)
Q_BINS = [0.0, 1e-4, 1e-3, 1e-2, 5e-2]


def main() -> None:
    assets = DEFAULT_RUNNER.run([DECODE_SMP_AP4M1_RAW_FILE, DECODE_VARIANT_INFO])
    raw, vi = assets[DECODE_SMP_AP4M1_RAW_FILE.asset_id], assets[DECODE_VARIANT_INFO.asset_id]
    assert isinstance(raw, FileAsset) and isinstance(vi, FileAsset)
    ea, oa = pl.col("effectAllele"), pl.col("otherAllele")
    snvs = (
        pl.read_csv(
            raw.path,
            separator="\t",
            columns=["Chrom", "Pos", "effectAllele", "otherAllele", "ImpMAF"],
            null_values=["."],
            quote_char=None,
            schema_overrides={"Pos": pl.Int64},
        )
        .with_columns(
            chrom_code_expr("Chrom").alias("CHR"),
            pl.len().over("Chrom", "Pos").alias("n_at_pos"),
        )
        .filter(
            (ea.str.len_bytes() == 1) & (oa.str.len_bytes() == 1) & (ea != oa) & (ea != "*")
        )
        .select(
            "CHR",
            pl.col("Pos").alias("POS"),
            ea.alias("EA"),
            "ImpMAF",
            (pl.col("n_at_pos") == 1).alias("lone"),
            (pl.col("ImpMAF") >= 0.01).fill_null(False).alias("common"),
        )
    )
    vi_at = (
        pl.scan_parquet(vi.path)
        .select("CHR", "POS", "A1", "ImpFreqA1")
        .unique()
        .join(snvs.select("CHR", "POS").unique().lazy(), on=["CHR", "POS"], how="semi")
        .collect()
    )
    q = (
        snvs.join(vi_at, on=["CHR", "POS"], how="left")
        .group_by("CHR", "POS", "EA", "lone", "common")
        .agg(
            pl.col("ImpFreqA1").filter(pl.col("A1") != pl.col("EA")).sum().alias("q"),
            pl.col("A1").filter(pl.col("A1") != pl.col("EA")).alias("other_alleles"),
            (pl.col("A1") == pl.col("EA")).any().alias("ea_in_vi"),
        )
    )
    print(f"SNV rows: {q.height}; EA missing from variant_info: {int((~q['ea_in_vi']).sum())}")
    q = q.with_columns(
        pl.col("q").cut(Q_BINS, left_closed=True, include_breaks=False).alias("q_bin")
    )
    summary = (
        q.group_by("lone", "common", "q_bin")
        .len()
        .with_columns((pl.col("len") / pl.col("len").sum().over("lone", "common")).alias("share"))
        .sort("lone", "common", "q_bin", descending=[True, True, False])
    )
    print(summary)
    print("\nlargest q among common lone SNVs:")
    print(
        q.filter("lone", "common")
        .sort("q", descending=True)
        .head(8)
        .select("CHR", "POS", "EA", "q", "other_alleles")
    )


if __name__ == "__main__":
    main()
