"""Independently validate the orientation changes the DecodeME rsID-assignment chain (gwaslab
harmonize + ANNOVAR join) made, using UKB-PPP (RABGAP1L, hg19-native IDs, UK European) allele
frequencies as an external truth for which allele is which.

For each DecodeME row that the chain changed (swapped EA/NEA, or strand-complemented), and a
random control sample of unchanged rows, look up the PPP row at the same hg19 position with the
same unordered allele pair, compute the PPP frequency of DecodeME's (post-chain) EA, and compare to
DecodeME's post-chain EAF. |diff| < 0.05 => orientation consistent; for palindromic SNPs the
swapped-vs-strand-flip interpretation is exactly what this discriminates (unless MAF ~ 0.5).
Also compare the PRE-chain EA/EAF, to see whether the change fixed or broke orientation.
"""
from pathlib import Path

import polars as pl

from mecfs_bio.build_system.task.ppp_database.allele_key import unordered_allele_key

S = Path("assets/base_asset_store")
D = S / "gwas/ME_CFS/DecodeME/processed"
PPP = S / "reference_data/ukbb_ppp/sun_et_al_2023/extracted/ukbb_ppp_rabgap1l_sumstats_stacked.parquet"
cols = ["SNPID", "CHR", "POS", "EA", "NEA", "EAF", "BETA"]


def load(name: str) -> pl.DataFrame:
    return pl.scan_parquet(D / name).select(cols).with_columns(pl.col(pl.Categorical).cast(pl.Utf8)).collect()


before = load("decode_me_gwas_1_liftover_to_37_parquet_file.parquet")
after = load("decode_me_gwas_1_assign_rsids_via_dbsnp150.parquet")
j = after.join(before.select("SNPID", *[pl.col(c).alias(f"b_{c}") for c in cols[3:]]), on="SNPID")
j = j.with_columns(
    pl.when((pl.col("EA") == pl.col("b_EA")) & (pl.col("NEA") == pl.col("b_NEA"))).then(pl.lit("unchanged"))
    .when((pl.col("EA") == pl.col("b_NEA")) & (pl.col("NEA") == pl.col("b_EA"))).then(pl.lit("swapped"))
    .otherwise(pl.lit("complemented/other")).alias("change"),
    ((pl.col("EA").str.len_chars() > 1) | (pl.col("NEA").str.len_chars() > 1)).alias("is_indel"),
    ((pl.col("BETA") - pl.col("b_BETA")).abs() < 1e-9).alias("beta_unchanged"),
    pl.col("EA").str.replace_many(["A", "C", "G", "T"], ["t", "g", "c", "a"]).str.to_uppercase().eq(pl.col("NEA")).alias("palindromic"),
)
changed = j.filter(pl.col("change") != "unchanged")
control = j.filter(pl.col("change") == "unchanged").sample(200_000, seed=0)
targets = pl.concat([changed, control]).with_columns(unordered_allele_key("EA", "NEA").alias("k"))
print(targets.group_by("change", "is_indel", "palindromic", "beta_unchanged").len().sort("change", "is_indel", "palindromic"))

ppp = (
    pl.scan_parquet(PPP)
    .select(
        pl.col("CHROM").cast(pl.Int64).alias("CHR"),
        pl.col("ID").str.split(":").list.get(1).cast(pl.Int64).alias("POS"),
        pl.col("ALLELE0").alias("p_A0"),
        pl.col("ALLELE1").alias("p_A1"),
        pl.col("A1FREQ").alias("p_A1FREQ"),
    )
    .with_columns(unordered_allele_key("p_A0", "p_A1").alias("k"))
    .join(targets.lazy().select("CHR", "POS", "k").unique(), on=["CHR", "POS", "k"], how="semi")
    .collect()
)
# Keep only unambiguous PPP matches (a single PPP row for the site/key).
ppp = ppp.filter(pl.len().over(["CHR", "POS", "k"]) == 1)
m = targets.join(ppp, on=["CHR", "POS", "k"], how="inner").with_columns(
    pl.when(pl.col("EA") == pl.col("p_A1")).then(pl.col("p_A1FREQ")).otherwise(1 - pl.col("p_A1FREQ")).alias("ppp_freq_of_EA"),
    pl.when(pl.col("b_EA") == pl.col("p_A1")).then(pl.col("p_A1FREQ"))
    .when(pl.col("b_EA") == pl.col("p_A0")).then(1 - pl.col("p_A1FREQ"))
    .otherwise(None).alias("ppp_freq_of_b_EA"),
).with_columns(
    ((pl.col("EAF") - pl.col("ppp_freq_of_EA")).abs() < 0.05).alias("after_consistent"),
    ((pl.col("b_EAF") - pl.col("ppp_freq_of_b_EA")).abs() < 0.05).alias("before_consistent"),
    ((pl.col("EAF") - 0.5).abs() < 0.08).alias("maf_near_half"),
)
with pl.Config(tbl_rows=40, tbl_cols=-1, tbl_width_chars=250):
    print(
        m.group_by("change", "is_indel", "palindromic", "beta_unchanged")
        .agg(
            pl.len().alias("n_matched_in_ppp"),
            pl.col("after_consistent").mean().alias("frac_after_consistent"),
            pl.col("before_consistent").mean().alias("frac_before_consistent"),
            pl.col("maf_near_half").mean().alias("frac_maf_near_half"),
        )
        .sort("change", "is_indel", "palindromic", "beta_unchanged")
    )
    bad = m.filter((pl.col("change") != "unchanged") & ~pl.col("after_consistent") & ~pl.col("maf_near_half"))
    print(f"changed rows inconsistent with PPP after the chain (excluding MAF~0.5): {bad.height}")
    print(bad.select("SNPID", "POS", "b_EA", "b_NEA", "b_EAF", "EA", "NEA", "EAF", "p_A0", "p_A1", "p_A1FREQ", "beta_unchanged").head(15))
