"""For key variants at a DecodeME locus, show the prior weight (relative to the locus
median) and PIP under uniform / Approach 1 / Approach 2, and each approach's top
per-family contrasts at those variants.

Run: pixi r python experiments/claude/polyfun_approach2_real_run/inspect_locus_variants.py \
    chr1_173500000_174500000_palindromes_keep 173855298 173878862 173838788 \
    | tee experiments/claude/polyfun_approach2_real_run/inspect_locus_variants_chr1.log
"""

import sys
from pathlib import Path

import polars as pl

BASE = Path("assets/base_asset_store/gwas/ME_CFS/DecodeME/analysis")
STEMS = {
    "approach1": "decode_me_polyfun_explain{locus}_l10",
    "approach2": "decode_me_polyfun_explain_l2_sldsc_prior{locus}_l10",
}
KEY = ["CHR", "POS", "EA", "NEA"]


def main(locus: str, positions: list[int]) -> None:
    uniform = pl.read_parquet(
        BASE / (STEMS["approach1"].format(locus=locus) + "_susie_uniform") / "combined_cs.parquet"
    ).select(*KEY, pl.col("PIP").alias("pip_uniform"))
    for name, stem in STEMS.items():
        run = BASE / (stem.format(locus=locus) + "_susie_polyfun")
        prior = pl.read_parquet(run / "prior.parquet")
        prior = prior.with_columns(
            (pl.col("prior_weight") / pl.col("prior_weight").median()).alias("prior_vs_median")
        )
        pips = pl.read_parquet(run / "combined_cs.parquet").select(*KEY, pl.col("PIP").alias("pip"))
        rows = (
            prior.filter(pl.col("POS").is_in(positions))
            .join(pips, on=KEY, how="left")
            .join(uniform, on=KEY, how="left")
            .select(*KEY, "prior_vs_median", "pip_uniform", "pip")
            .sort("POS")
        )
        print(f"\n===== {name}: prior (x locus median) and PIP =====")
        print(rows)
        fam = pl.read_parquet(BASE / (stem.format(locus=locus) + "_explain_contrast") / "per_family_contrast.parquet")
        pos_col = "POS" if "POS" in fam.columns else fam.columns[1]
        print(f"--- {name}: top family contrasts at these variants ---")
        for pos in positions:
            sub = fam.filter(pl.col(pos_col) == pos)
            val = [c for c in sub.columns if c not in (*KEY, "family")][0]
            print(pos, sub.sort(val, descending=True).select("family", val).head(3).rows())


if __name__ == "__main__":
    main(sys.argv[1], [int(p) for p in sys.argv[2:]])
