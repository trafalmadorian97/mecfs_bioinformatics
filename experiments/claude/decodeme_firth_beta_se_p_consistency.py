"""
Check whether BETA/SE is consistent with CHISQ and LOG10P in the DecodeME GWAS 1 QC-filtered summary statistics.

REGENIE's docs warn that with Firth correction the reported SE may not be consistent with the p-value.
We compare three quantities per variant:
  z2_wald = (BETA/SE)^2
  CHISQ   (reported test statistic)
  chi2 implied by LOG10P (1-df)

Run: pixi r python experiments/claude/decodeme_firth_beta_se_p_consistency.py | tee experiments/claude/decodeme_firth_beta_se_p_consistency.log
"""

import numpy as np
import polars as pl
from scipy import stats

PATH = "assets/base_asset_store/gwas/ME_CFS/DecodeME/processed/decode_me_gwas_1_filtered_for_quality_control.parquet"

pl.Config.set_tbl_cols(20)
pl.Config.set_tbl_width_chars(250)

df = pl.read_parquet(
    PATH,
    columns=[
        "ID",
        "A1FREQ",
        "N_CASES",
        "TEST",
        "BETA",
        "SE",
        "CHISQ",
        "LOG10P",
        "EXTRA",
    ],
)
print("rows:", df.height)
print("\nTEST values:\n", df["TEST"].value_counts())
print("\nEXTRA values (top):\n", df["EXTRA"].value_counts(sort=True).head(10))

log10p = df["LOG10P"].to_numpy()
# chi2 (1 df) implied by LOG10P; isf of p = 10^-log10p, computed in log space for tiny p
logp = -log10p * np.log(10)
chi2_from_p = stats.chi2.isf(np.exp(logp), df=1)
# for very small p use the log-space inverse
small = log10p > 250
if small.any():
    chi2_from_p[small] = [
        stats.norm.isf(np.exp(lp - np.log(2))) ** 2 for lp in logp[small]
    ]

df = df.with_columns(
    z2_wald=(pl.col("BETA") / pl.col("SE")) ** 2,
    chi2_from_p=pl.Series(chi2_from_p),
).with_columns(
    ratio_wald_vs_chisq=pl.col("z2_wald") / pl.col("CHISQ"),
    ratio_chisq_vs_p=pl.col("CHISQ") / pl.col("chi2_from_p"),
    ratio_wald_vs_p=pl.col("z2_wald") / pl.col("chi2_from_p"),
    is_firth=pl.col("EXTRA").str.contains("(?i)firth")
    | pl.col("EXTRA").str.contains("SE_"),
)


def summarize(sub: pl.DataFrame, label: str) -> None:
    print(f"\n=== {label}: n={sub.height} ===")
    if sub.height == 0:
        return
    for col in ["ratio_wald_vs_chisq", "ratio_chisq_vs_p", "ratio_wald_vs_p"]:
        s = sub[col]
        rel = (s - 1).abs()
        print(
            f"{col:22s} median={s.median():.6f} min={s.min():.6f} max={s.max():.6f} "
            f"frac|r-1|>1e-3={(rel > 1e-3).mean():.5f} frac|r-1|>0.01={(rel > 0.01).mean():.5f} "
            f"frac|r-1|>0.1={(rel > 0.1).mean():.5f}"
        )


summarize(df, "all variants")
summarize(df.filter(pl.col("LOG10P") > 5), "LOG10P > 5")
summarize(df.filter(pl.col("CHISQ") > 1), "CHISQ > 1")

# stratify by EXTRA flag
for extra_val in df["EXTRA"].unique().to_list()[:10]:
    summarize(df.filter(pl.col("EXTRA") == extra_val), f"EXTRA == {extra_val!r}")

print("\nWorst Wald-vs-p discrepancies among CHISQ>1:")
print(
    df.filter(pl.col("CHISQ") > 1)
    .with_columns(absdev=(pl.col("ratio_wald_vs_p") - 1).abs())
    .sort("absdev", descending=True)
    .head(15)
)

print("\nTop hits by LOG10P:")
print(df.sort("LOG10P", descending=True).head(15))

# discrepancy vs MAF
df = df.with_columns(
    maf=pl.min_horizontal(pl.col("A1FREQ"), 1 - pl.col("A1FREQ"))
).with_columns(
    maf_bin=pl.col("maf").cut([0.005, 0.01, 0.05, 0.1, 0.25]),
)
print("\nWald-vs-p discrepancy by MAF bin (CHISQ>1):")
print(
    df.filter(pl.col("CHISQ") > 1)
    .group_by("maf_bin")
    .agg(
        n=pl.len(),
        median_ratio=pl.col("ratio_wald_vs_p").median(),
        frac_dev_gt_1pct=((pl.col("ratio_wald_vs_p") - 1).abs() > 0.01).mean(),
        max_abs_dev=(pl.col("ratio_wald_vs_p") - 1).abs().max(),
    )
    .sort("maf_bin")
)
