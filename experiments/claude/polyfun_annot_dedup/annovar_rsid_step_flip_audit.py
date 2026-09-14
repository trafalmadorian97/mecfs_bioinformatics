"""Genome-wide: how often does the DecodeME ANNOVAR rsID-assignment step swap EA/NEA relative to its
input (the build-37 liftover table)? For swapped rows, check BETA/EAF/cases/controls consistency and
which orientation matches the hg19 FASTA; also check its REF/ALT columns against the FASTA."""
from pathlib import Path

import numpy as np
import polars as pl

from experiments.claude.polyfun_annot_dedup.benchmark_fasta_ref_check import HG19, read_fai, ref_match

D = Path("assets/base_asset_store/gwas/ME_CFS/DecodeME/processed")
cols = ["SNPID", "CHR", "POS", "EA", "NEA", "EAF", "BETA", "A1FREQ_CASES", "A1FREQ_CONTROLS"]


def load(name: str, extra: list[str]) -> pl.DataFrame:
    lf = pl.scan_parquet(D / name)
    have = lf.collect_schema().names()
    return lf.select([c for c in cols + extra if c in have]).with_columns(pl.col(pl.Categorical).cast(pl.Utf8)).collect()


before = load("decode_me_gwas_1_liftover_to_37_parquet_file.parquet", [])
after = load("decode_me_gwas_1_assign_rsids_via_dbsnp150.parquet", ["REF", "ALT", "rsID"])
print(f"liftover-37 rows={before.height:,} dup SNPID={before.height - before['SNPID'].n_unique():,}")
print(f"rsid-assigned rows={after.height:,} dup SNPID={after.height - after['SNPID'].n_unique():,}")
print("rsid-assigned columns:", after.columns)

j = after.join(before.select("SNPID", *[pl.col(c).alias(f"b_{c}") for c in cols if c != "SNPID"]), on="SNPID", how="inner")
same = (j["EA"] == j["b_EA"]) & (j["NEA"] == j["b_NEA"])
swap = (j["EA"] == j["b_NEA"]) & (j["NEA"] == j["b_EA"])
indel = (j["b_EA"].str.len_chars() > 1) | (j["b_NEA"].str.len_chars() > 1)
print(f"\njoined rows={j.height:,}  same orientation={int(same.sum()):,}  swapped={int(swap.sum()):,}  other={int((~same & ~swap).sum()):,}")
print(f"swapped SNVs={int((swap & ~indel).sum()):,}  swapped indels={int((swap & indel).sum()):,}  (indels total={int(indel.sum()):,})")

s = j.filter(swap)
if s.height:
    print("swapped: BETA negated:", bool(np.allclose(s["BETA"], -s["b_BETA"])))
    print("swapped: EAF complemented:", bool(np.allclose(s["EAF"], 1 - s["b_EAF"], atol=1e-5)))
    for c in ["A1FREQ_CASES", "A1FREQ_CONTROLS"]:
        if c in s.columns:
            print(f"swapped: {c} complemented: {bool(np.allclose(s[c], 1 - s[f'b_{c}'], atol=1e-5))}, unchanged: {bool(np.allclose(s[c], s[f'b_{c}'], atol=1e-5))}")
    fai = read_fai(HG19)
    mm = np.memmap(HG19, dtype=np.uint8, mode="r")
    ch, ps = s["CHR"].cast(pl.Utf8).to_numpy(), s["POS"].to_numpy()
    b_nea = ref_match(mm, fai, ch, ps, s["b_NEA"])
    a_nea = ref_match(mm, fai, ch, ps, s["NEA"])
    longer_in = np.where(s["b_NEA"].str.len_chars() >= s["b_EA"].str.len_chars(), "NEA longer (deletion in input)", "EA longer (insertion in input)")
    print(f"swapped: input NEA==hg19 ref {int(b_nea.sum()):,}/{s.height:,}; output NEA==hg19 ref {int(a_nea.sum()):,}/{s.height:,}")
    print(pl.Series("input indel type", longer_in).value_counts())
    with pl.Config(tbl_cols=-1, tbl_width_chars=250):
        print(s.select("SNPID", "POS", "b_EA", "b_NEA", "EA", "NEA", *[c for c in ["REF", "ALT", "rsID"] if c in s.columns]).head(10))

if "REF" in after.columns:
    fai = read_fai(HG19)
    mm = np.memmap(HG19, dtype=np.uint8, mode="r")
    a = after.filter(pl.col("REF").is_not_null())
    m = ref_match(mm, fai, a["CHR"].cast(pl.Utf8).to_numpy(), a["POS"].to_numpy(), a["REF"])
    ind = ((a["REF"].str.len_chars() > 1) | (a["ALT"].str.len_chars() > 1)).to_numpy()
    print(f"\nREF column == hg19 ref: SNVs {m[~ind].mean():.4%} (n={int((~ind).sum()):,}), indels {m[ind].mean():.4%} (n={int(ind.sum()):,})")
    eq = (a["NEA"] == a["REF"]).to_numpy()
    print(f"NEA == REF column: SNVs {eq[~ind].mean():.4%}, indels {eq[ind].mean():.4%}")

print("\n=== Per-row breakdown of swapped rows ===")
tol = 1e-4
s = s.with_columns(
    pl.when(pl.col("BETA").is_nan() | pl.col("b_BETA").is_nan() | pl.col("BETA").is_null()).then(pl.lit("nan/null"))
    .when((pl.col("BETA") + pl.col("b_BETA")).abs() <= tol * (1 + pl.col("b_BETA").abs())).then(pl.lit("negated"))
    .when((pl.col("BETA") - pl.col("b_BETA")).abs() <= tol * (1 + pl.col("b_BETA").abs())).then(pl.lit("unchanged"))
    .otherwise(pl.lit("other")).alias("beta_status"),
    pl.when(pl.col("EAF").is_nan() | pl.col("b_EAF").is_nan() | pl.col("EAF").is_null()).then(pl.lit("nan/null"))
    .when((pl.col("EAF") - (1 - pl.col("b_EAF"))).abs() <= 1e-4).then(pl.lit("complemented"))
    .when((pl.col("EAF") - pl.col("b_EAF")).abs() <= 1e-4).then(pl.lit("unchanged"))
    .otherwise(pl.lit("other")).alias("eaf_status"),
    ((pl.col("b_EA").str.len_chars() > 1) | (pl.col("b_NEA").str.len_chars() > 1)).alias("is_indel"),
    pl.Series("input_nea_is_ref", b_nea),
)
with pl.Config(tbl_rows=40, tbl_cols=-1, tbl_width_chars=250):
    print(s.group_by("is_indel", "input_nea_is_ref", "beta_status", "eaf_status").len().sort("is_indel", "input_nea_is_ref", "beta_status"))
    print(s.filter(pl.col("beta_status") == "other").select("SNPID", "b_EA", "b_NEA", "EA", "NEA", "b_BETA", "BETA", "b_EAF", "EAF").head(5))

print("\n=== 'other' rows (neither same nor swapped) ===")
o = j.filter(~same & ~swap)
with pl.Config(tbl_cols=-1, tbl_width_chars=250):
    print(o.select("SNPID", "POS", "b_POS", "b_EA", "b_NEA", "EA", "NEA", "REF", "ALT", "b_BETA", "BETA").head(10))
    print("POS changed:", int((o["POS"] != o["b_POS"]).sum()), "of", o.height)

print("\n=== Are the beta-unchanged swaps and the 'other' rows strand flips? ===")
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def revcomp(expr: pl.Expr) -> pl.Expr:
    return expr.str.reverse().str.replace_many(["A", "C", "G", "T"], ["t", "g", "c", "a"]).str.to_uppercase()


pal = s.filter(~pl.col("is_indel")).with_columns(
    (revcomp(pl.col("b_EA")) == pl.col("b_NEA")).alias("palindromic")
)
with pl.Config(tbl_rows=20):
    print(pal.group_by("beta_status", "palindromic").len().sort("beta_status"))
o2 = o.with_columns(
    ((revcomp(pl.col("b_EA")) == pl.col("EA")) & (revcomp(pl.col("b_NEA")) == pl.col("NEA"))).alias("strand_complement"),
    (pl.col("BETA") == pl.col("b_BETA")).alias("beta_unchanged"),
)
print(o2.group_by("strand_complement", "beta_unchanged").len())
unch = pal.filter(pl.col("beta_status") == "unchanged")
print("beta-unchanged swaps, chromosome/position spread (hg19):")
print(unch.group_by("CHR").agg(pl.len(), pl.col("POS").min(), pl.col("POS").max()).sort("CHR").head(25))
