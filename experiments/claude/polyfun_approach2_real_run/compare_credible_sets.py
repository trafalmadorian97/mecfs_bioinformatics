"""Compare SUSIE credible sets at a DecodeME locus between the precomputed PolyFun
prior (Approach 1) and the DecodeME L2-regularized S-LDSC prior (Approach 2), plus
the uniform-prior baseline, for each run config.

Run: pixi r python experiments/claude/polyfun_approach2_real_run/compare_credible_sets.py \
    chr1_173500000_174500000_palindromes_keep \
    | tee experiments/claude/polyfun_approach2_real_run/compare_credible_sets_chr1.log
"""

import sys
from pathlib import Path

import polars as pl

BASE = Path("assets/base_asset_store/gwas/ME_CFS/DecodeME/analysis")
RUNS = {
    "uniform": "decode_me_polyfun_explain{locus}_{cfg}_susie_uniform",
    "approach1": "decode_me_polyfun_explain{locus}_{cfg}_susie_polyfun",
    "approach2": "decode_me_polyfun_explain_l2_sldsc_prior{locus}_{cfg}_susie_polyfun",
}
CONFIGS = ("l1", "l2", "l10", "l10_strict")
KEY = ["CHR", "POS", "EA", "NEA"]


def _cs(run: str, cfg: str, locus: str) -> pl.DataFrame:
    path = BASE / RUNS[run].format(locus=locus, cfg=cfg) / "combined_cs.parquet"
    return pl.read_parquet(path)


def _describe(cs: pl.DataFrame) -> str:
    if cs.height == 0:
        return "no credible sets"
    parts = []
    for name, grp in cs.group_by("cs", maintain_order=True):
        top = grp.sort("PIP", descending=True).row(0, named=True)
        parts.append(
            f"{name[0]}: {grp.height} vars, top {top['CHR']}:{top['POS']}:{top['NEA']}:{top['EA']} PIP={top['PIP']:.3f}"
        )
    return "; ".join(parts)


def _variants(cs: pl.DataFrame) -> set[tuple]:
    return set(cs.select(KEY).iter_rows())


def main(locus: str) -> None:
    for cfg in CONFIGS:
        print(f"\n===== {cfg} =====")
        sets = {run: _cs(run, cfg, locus) for run in RUNS}
        for run, cs in sets.items():
            print(f"  {run:>9}: {_describe(cs)}")
        v = {run: _variants(cs) for run, cs in sets.items()}
        for a, b in (("approach1", "approach2"), ("uniform", "approach1"), ("uniform", "approach2")):
            inter, union = v[a] & v[b], v[a] | v[b]
            jac = len(inter) / len(union) if union else float("nan")
            print(f"  CS variants {a} vs {b}: |A|={len(v[a])} |B|={len(v[b])} shared={len(inter)} jaccard={jac:.2f}")
        pips = (
            sets["approach1"].select(*KEY, pl.col("PIP").alias("pip_a1"))
            .join(sets["approach2"].select(*KEY, pl.col("PIP").alias("pip_a2")), on=KEY, how="full", coalesce=True)
            .fill_null(0.0)
            .with_columns((pl.col("pip_a2") - pl.col("pip_a1")).alias("delta"))
            .sort(pl.col("delta").abs(), descending=True)
        )
        print("  largest CS-variant PIP changes (approach2 - approach1):")
        print(pips.head(5))


if __name__ == "__main__":
    main(sys.argv[1])
