"""Inspect the mirrored-SNV collisions in the RABGAP1L template protein and whether they are in the
PPP HapMap3 index (where align_protein_to_index's arbitrary unique() could pick the wrong row)."""
from pathlib import Path

import polars as pl

from mecfs_bio.build_system.task.ppp_database.allele_key import unordered_allele_key

STORE = Path("assets/base_asset_store")
T = STORE / "reference_data/ukbb_ppp/sun_et_al_2023/extracted/ukbb_ppp_rabgap1l_sumstats_stacked.parquet"
IDX = STORE / "reference_data/ukbb_ppp_variant_index/hapmap_3_membership_list/processed/ppp_variant_index.parquet"
for chrom in (7, 11):
    df = (
        pl.scan_parquet(T)
        .filter(pl.col("CHROM") == chrom)
        .filter((pl.col("ALLELE0").str.len_chars() == 1) & (pl.col("ALLELE1").str.len_chars() == 1))
        .select("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "BETA", "SE")
        .with_columns(unordered_allele_key("ALLELE0", "ALLELE1").alias("k"))
        .filter(pl.len().over(["GENPOS", "k"]) > 1)
        .collect()
    )
    with pl.Config(tbl_cols=-1, tbl_width_chars=250):
        print(df)
        for pos in df["GENPOS"].unique():
            print("in PPP index:", pl.scan_parquet(IDX).filter((pl.col("CHR") == chrom) & (pl.col("POS") == pos)).collect())
