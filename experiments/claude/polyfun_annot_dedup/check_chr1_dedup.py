"""Run the current BuildBaselineLFAnnotationParquetTask dedup on the locally cached chr1
member and characterize any (CHR, BP, unordered-allele-key) groups with conflicting annotations."""
from pathlib import Path

import polars as pl

from mecfs_bio.build_system.task.annotation_weights.build_baseline_lf_annotation_parquet_task import (
    ANNOT_KEY_COLUMNS,
    _dedup_one_chromosome,
)
from mecfs_bio.build_system.task.ppp_database.allele_key import unordered_allele_key

MEMBER = Path(
    "assets/base_asset_store/reference_data/polyfun/annotations/raw/"
    "baseline_lf_2.2_ukb_annot_parquet_members/baselineLF2.2.UKB.1.annot.parquet"
)
try:
    _dedup_one_chromosome(MEMBER)
    print("dedup PASSED on local member")
except AssertionError as e:
    print("dedup FAILED:", e)

lazy = pl.scan_parquet(MEMBER)
annot = [c for c in lazy.collect_schema().names() if c not in ANNOT_KEY_COLUMNS]
df = lazy.with_columns(unordered_allele_key("A1", "A2").alias("k")).collect()
dup = df.filter(pl.len().over(["CHR", "BP", "k"]) > 1).sort("BP")
print("rows in dup key groups:", dup.height)
# Raw (float64) vs float32 conflicts
for label, frame in [("raw", dup), ("f32", dup.with_columns([pl.col(c).cast(pl.Float32) for c in annot]))]:
    u = frame.unique(subset=["CHR", "BP", "k", *annot])
    print(label, "conflicting groups:", u.height - u.n_unique(subset=["CHR", "BP", "k"]))
print(dup.select(["CHR", "BP", "SNP", "A1", "A2"]).head(20))
# which columns differ in conflicting groups
conf = dup.unique(subset=["CHR", "BP", "k", *annot])
conf = conf.filter(pl.len().over(["CHR", "BP", "k"]) > 1)
diff_cols = [c for c in annot if conf.group_by(["BP", "k"]).agg(pl.col(c).n_unique().alias("n")).filter(pl.col("n") > 1).height > 0]
print("differing columns:", diff_cols)
print(conf.select(["BP", "SNP", "A1", "A2", *diff_cols[:6]]).sort("BP").head(20))
print("dtypes:", {c: str(df.schema[c]) for c in annot[:3]})
nan_cols = [c for c in annot if df[c].is_nan().any()] if df.schema[annot[0]].is_float() else []
print("columns with NaN:", nan_cols[:10])

# Per-pair detail: how many annotation columns differ, and are all conflicting pairs indels?
pairs = dup.group_by(["BP", "k"]).agg(
    [pl.col("A1").str.len_chars().max().alias("maxlen")]
    + [(pl.col(c).n_unique() > 1).alias(c) for c in annot]
)
pairs = pairs.with_columns(pl.sum_horizontal(annot).alias("n_diff"))
print("all conflicting pairs are indels:", (pairs["maxlen"] > 1).all())
print("n_diff distribution:", pairs["n_diff"].describe())
ex = dup.filter(pl.col("BP") == 4372921)
print(ex.transpose(include_header=True).filter(pl.col("column_0") != pl.col("column_1")))
# Indel swapped-allele pairs that DID agree (collapsed silently)
all_pairs = df.filter(pl.len().over(["CHR", "BP", "k"]) > 1)
print("total rows in swapped groups (incl agreeing):", all_pairs.height)
