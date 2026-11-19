"""How many candidate (locus, protein) pairs were dropped by the selection budget in the pain/fibro PPP coloc sweeps."""

import polars as pl

pl.Config.set_tbl_rows(60)
pl.Config.set_tbl_width_chars(200)

RUNS = {
    "fibro": "assets/base_asset_store/gwas/fibromyalgia/kerrebijn_et_al/analysis/kerrebijn_fibro_ppp_coloc_full_ppp_coloc/raw_results.parquet",
    "pain": "assets/base_asset_store/gwas/multisite_pain/johnston_et_al/analysis/johnston_et_al_pain_ppp_coloc_full_ppp_coloc/raw_results.parquet",
}

for name, path in RUNS.items():
    d = pl.read_parquet(path)
    print(f"\n===== {name} =====")
    print(d["STATUS"].value_counts().sort("count", descending=True))
    pairs = d.filter(pl.col("PROTEIN_ID").is_not_null())
    print("candidate pairs:", pairs.height)
    per_locus = (
        pairs.group_by("LOCUS_ID")
        .agg(
            pl.len().alias("n_cand"),
            (pl.col("STATUS") == "not_selected_budget").sum().alias("n_dropped"),
            (pl.col("STATUS") == "tested").sum().alias("n_tested"),
            (pl.col("PP_H4") >= 0.8).sum().alias("n_h4_ge_0.8"),
            pl.col("PP_H4").max().alias("max_h4"),
            (pl.col("CIS_STATUS") == "cis").sum().alias("n_cis"),
            ((pl.col("CIS_STATUS") == "cis") & (pl.col("STATUS") == "not_selected_budget")).sum().alias("n_cis_dropped"),
            pl.col("LEAD_MLOG10P").first(),
        )
        .sort("n_dropped", descending=True)
    )
    print("loci:", per_locus.height, " loci with drops:", (per_locus["n_dropped"] > 0).sum())
    print(per_locus.filter(pl.col("n_dropped") > 0))
    print("total dropped:", per_locus["n_dropped"].sum(), " cis dropped:", per_locus["n_cis_dropped"].sum())

print("\n===== loci with drops AND a tested PP.H4 >= 0.8; plus dropped cis pairs =====")
for name, path in RUNS.items():
    d = pl.read_parquet(path).filter(pl.col("PROTEIN_ID").is_not_null())
    hot = (
        d.group_by("LOCUS_ID")
        .agg(
            (pl.col("STATUS") == "not_selected_budget").sum().alias("n_dropped"),
            (pl.col("PP_H4") >= 0.8).sum().alias("n_tested_h4_ge_0.8"),
            pl.col("PP_H4").max().alias("max_h4"),
        )
        .filter((pl.col("n_dropped") > 0) & (pl.col("max_h4") >= 0.8))
        .sort("n_dropped", descending=True)
    )
    print(f"\n--- {name} ---")
    print(hot)
    print(
        d.filter((pl.col("STATUS") == "not_selected_budget") & (pl.col("CIS_STATUS") == "cis")).select(
            "LOCUS_ID", "GENE", "PROTEIN_MLOG10P", "SCORE"
        )
    )
