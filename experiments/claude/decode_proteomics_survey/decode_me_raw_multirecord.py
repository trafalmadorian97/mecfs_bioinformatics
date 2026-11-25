"""
Multi-record positions in DecodeME GWAS-1's raw REGENIE output, before and after DecodeME's own
QC variant list (gwas_qced.var.gz).

  pixi r python experiments/claude/decode_proteomics_survey/decode_me_raw_multirecord.py
"""

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.filtered_snps_gwas_1 import (
    DECODE_ME_FILTER_SNPS_GWAS_1_TASK,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.raw_gwas_data.decode_me_gwas_1 import (
    DECODE_ME_GWAS_1_TASK,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset

pl.Config.set_tbl_rows(16)
pl.Config.set_tbl_width_chars(200)


def summarize(label: str, frame: pl.DataFrame, chrom: str, pos: str, a0: str, a1: str) -> None:
    frame = frame.with_columns(pl.len().over(chrom, pos).alias("n_at_pos"))
    multi = frame.filter(pl.col("n_at_pos") > 1)
    positions = frame.select(pl.struct(chrom, pos).n_unique()).item()
    print(
        f"\n{label}: {frame.height} records at {positions} positions; {multi.height} records "
        f"({multi.height / frame.height:.2%}) at "
        f"{multi.select(pl.struct(chrom, pos).n_unique()).item()} multi-record positions"
    )
    print(multi.sort(chrom, pos).head(12).select(chrom, pos, "ID", a0, a1, "A1FREQ", "BETA", "SE"))


def main() -> None:
    for task in (DECODE_ME_GWAS_1_TASK, DECODE_ME_FILTER_SNPS_GWAS_1_TASK):
        asset = DEFAULT_RUNNER.run([task])[task.asset_id]
        frame = scan_dataframe_asset(asset, meta=task.meta, parquet_backend="polars").to_native()
        frame = frame.collect() if isinstance(frame, pl.LazyFrame) else frame
        print(f"{task.asset_id} columns: {frame.columns}")
        summarize(task.asset_id, frame, "CHROM", "GENPOS", "ALLELE0", "ALLELE1")


if __name__ == "__main__":
    main()
