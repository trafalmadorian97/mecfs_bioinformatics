"""
How many positions of the DecodeME GWAS-1 (GRCh38-harmonized) carry more than one record, i.e.
multi-allelic sites split into several rows?

  pixi r python experiments/claude/decode_proteomics_survey/trait_multirecord_positions.py
"""

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_build_38_harmonized import (
    DECODE_ME_GWAS_1_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset

task = DECODE_ME_GWAS_1_BUILD_38_HARMONIZED
asset = DEFAULT_RUNNER.run([task])[task.asset_id]
frame = (
    scan_dataframe_asset(asset, meta=task.meta, parquet_backend="polars")
    .to_native()
    .select("CHR", "POS", "EA", "NEA", "EAF")
    .collect()
)
ea, nea = pl.col("EA"), pl.col("NEA")
frame = frame.with_columns(
    pl.len().over("CHR", "POS").alias("n_at_pos"),
    ((ea.str.len_bytes() == 1) & (nea.str.len_bytes() == 1)).alias("snv"),
    (pl.min_horizontal("EAF", 1 - pl.col("EAF")) >= 0.01).fill_null(False).alias("common"),
)
positions = frame.select(pl.struct("CHR", "POS").n_unique()).item()
multi = frame.filter(pl.col("n_at_pos") > 1)
print(f"{frame.height} records at {positions} positions")
print(
    f"records at multi-record positions: {multi.height} ({multi.height / frame.height:.2%}); "
    f"positions: {multi.select(pl.struct('CHR', 'POS').n_unique()).item()}"
)
print(
    frame.group_by("snv", "common")
    .agg(pl.len(), (pl.col("n_at_pos") > 1).sum().alias("at_multi"), (pl.col("n_at_pos") > 1).mean().alias("share"))
    .sort("snv", "common")
)
pl.Config.set_tbl_rows(12)
print(multi.sort("CHR", "POS").head(12))
