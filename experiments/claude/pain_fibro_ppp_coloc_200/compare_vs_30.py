"""
Compare the 200-per-locus pain/fibro PPP coloc sweeps with the earlier 30-per-locus run. The
30-per-locus outputs were overwritten, so pairs that cap 30 would have dropped are approximated as
those ranked beyond 30 by SCORE within their locus (the selector's tie-breakers are ignored).
"""

import polars as pl

pl.Config.set_tbl_rows(60)
pl.Config.set_tbl_width_chars(220)
pl.Config.set_fmt_str_lengths(80)

RUNS = {
    "fibro": "assets/base_asset_store/gwas/fibromyalgia/kerrebijn_et_al/analysis/kerrebijn_fibro_ppp_coloc_full_ppp_coloc",
    "pain": "assets/base_asset_store/gwas/multisite_pain/johnston_et_al/analysis/johnston_et_al_pain_ppp_coloc_full_ppp_coloc",
}
OLD_CAP = 30

for name, base in RUNS.items():
    raw = pl.read_parquet(f"{base}/raw_results.parquet").filter(pl.col("PROTEIN_ID").is_not_null())
    display = pl.read_parquet(f"{base}_display.parquet")
    print(f"\n===== {name} =====")
    print(raw["STATUS"].value_counts().sort("count", descending=True))
    print("display rows:", display.height, " fdr_sig:", display["fdr_sig"].sum())
    ranked = raw.with_columns(
        pl.col("SCORE").rank("ordinal", descending=True).over("LOCUS_ID").alias("RANK_IN_LOCUS")
    )
    new = ranked.filter((pl.col("RANK_IN_LOCUS") > OLD_CAP) & (pl.col("STATUS") == "tested"))
    print("newly tested pairs (rank >30 in locus):", new.height, " with PP.H4>=0.8:", (new["PP_H4"] >= 0.8).sum())
    print(
        new.group_by("LOCUS_ID")
        .agg(
            pl.len().alias("n_new"),
            (pl.col("PP_H4") >= 0.8).sum().alias("n_new_h4_ge_0.8"),
            pl.col("PP_H4").max().alias("max_new_h4"),
            pl.col("GENE").filter(pl.col("PP_H4") >= 0.8).head(8).str.join(",").alias("new_hits"),
        )
        .sort("n_new_h4_ge_0.8", descending=True)
        .filter(pl.col("n_new_h4_ge_0.8") > 0)
    )
    print("cis among new hits:", new.filter((pl.col("PP_H4") >= 0.8) & (pl.col("CIS_STATUS") == "cis"))["GENE"].to_list())
    still = raw.filter(pl.col("STATUS") == "not_selected_budget")
    print("still dropped:", still.height, still.group_by("LOCUS_ID").len().sort("len", descending=True))
