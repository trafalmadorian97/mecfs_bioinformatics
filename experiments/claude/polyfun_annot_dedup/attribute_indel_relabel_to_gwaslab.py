"""Attribute the insertion<->deletion relabelling to a pipeline step: compare the liftover table,
the gwaslab harmonized dump, and the ANNOVAR-joined table for two example indels, and count
genome-wide how many indel swaps are already present in the gwaslab harmonized dump."""
from pathlib import Path

import polars as pl

D = Path("assets/base_asset_store/gwas/ME_CFS/DecodeME/processed")
tables = {
    "liftover_37": "decode_me_gwas_1_liftover_to_37_parquet_file.parquet",
    "gwaslab_harmonized_dump": "decode_me_gwas_1_harmonized_dump_to_parquet.parquet",
    "annovar_joined": "decode_me_gwas_1_assign_rsids_via_dbsnp150.parquet",
}
frames = {
    k: pl.scan_parquet(D / v).select("SNPID", "POS", "EA", "NEA", "EAF", "BETA", "STATUS" if "STATUS" in pl.scan_parquet(D / v).collect_schema().names() else pl.lit(None).alias("STATUS"))
    .with_columns(pl.col(pl.Categorical).cast(pl.Utf8)).collect()
    for k, v in tables.items()
}
with pl.Config(tbl_cols=-1, tbl_width_chars=200):
    for snpid in ["15:54405994:TG:T", "1:7857852:G:GA"]:
        for k, f in frames.items():
            print(k, f.filter(pl.col("SNPID") == snpid).rows())
b, h = frames["liftover_37"], frames["gwaslab_harmonized_dump"]
j = h.join(b.select("SNPID", pl.col("EA").alias("b_EA"), pl.col("NEA").alias("b_NEA")), on="SNPID")
indel = (pl.col("b_EA").str.len_chars() > 1) | (pl.col("b_NEA").str.len_chars() > 1)
swapped = (pl.col("EA") == pl.col("b_NEA")) & (pl.col("NEA") == pl.col("b_EA"))
print(j.select(swapped.and_(indel).sum().alias("indel_swaps_in_gwaslab_dump"), swapped.and_(indel.not_()).sum().alias("snv_swaps_in_gwaslab_dump"), pl.len()))
