"""Count, per chromosome, (CHR, BP, unordered-allele-key) groups with >1 row in the baseline-LF
annotation members and in the PolyFun precomputed snpvar prior; and whether any fall in the
DecodeME chr1:173.5M-174.5M test locus."""
from pathlib import Path

import polars as pl

from mecfs_bio.build_system.task.annotation_weights.build_baseline_lf_annotation_parquet_task import (
    ANNOT_KEY_COLUMNS,
)
from mecfs_bio.build_system.task.ppp_database.allele_key import unordered_allele_key

MEMBERS = Path("assets/base_asset_store/reference_data/polyfun/annotations/raw/baseline_lf_2.2_ukb_annot_parquet_members")
PRIOR = Path("assets/base_asset_store/reference_data/polyfun/precomputed_prior/raw/polyfun_precomputed_heritability_weight_concat.parquet")
K = ["CHR", "BP", "k"]

tot_groups = tot_conf = 0
for chrom in range(1, 23):
    lazy = pl.scan_parquet(MEMBERS / f"baselineLF2.2.UKB.{chrom}.annot.parquet")
    annot = [c for c in lazy.collect_schema().names() if c not in ANNOT_KEY_COLUMNS]
    dup = (
        lazy.with_columns(unordered_allele_key("A1", "A2").alias("k"))
        .filter(pl.len().over(K) > 1)
        .collect()
    )
    groups = dup.n_unique(subset=K)
    u = dup.unique(subset=[*K, *annot])
    conf = u.height - u.n_unique(subset=K)
    tot_groups += groups
    tot_conf += conf
    in_locus = dup.filter((pl.col("CHR") == 1) & pl.col("BP").is_between(173_500_000, 174_500_000)).height
    print(f"chr{chrom}: swapped groups={groups} conflicting={conf} rows_in_test_locus={in_locus}")
print("TOTAL annot swapped groups", tot_groups, "conflicting", tot_conf)

p = pl.scan_parquet(PRIOR)
print(p.collect_schema())
pdup = (
    p.with_columns(unordered_allele_key("A1", "A2").alias("k"))
    .filter(pl.len().over(K) > 1)
    .collect()
)
print("prior rows in swapped groups:", pdup.height, "groups:", pdup.n_unique(subset=K))
print(pdup.sort("CHR", "BP").head(6))
print("prior swapped rows in test locus:", pdup.filter((pl.col("CHR") == 1) & pl.col("BP").is_between(173_500_000, 174_500_000)).height)
