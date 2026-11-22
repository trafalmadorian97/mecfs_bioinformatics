"""
Can ImpMAF pick the right variant_info row when the effect allele alone does not?

For every row of one protein file, candidates are the deduplicated variant_info rows at the same
(Chrom, Pos) with A1 == effectAllele. A candidate agrees when |min(f, 1 - f) - ImpMAF| <= 1e-5,
with f = ImpFreqA1. The rule under test: EAF = f when the agreeing candidates all have the same f
(within 1e-5), else EAF is unknown (null). The script counts each case.

Reads the parquet caches written by check_variant_info_alleles.py.

Usage (from the repo root):
  pixi r python experiments/claude/decode_proteomics_survey/check_variant_info_maf_disambiguation.py \
      <cache_dir> | tee experiments/claude/decode_proteomics_survey/check_variant_info_maf_disambiguation.log
"""

import sys
from pathlib import Path

import polars as pl

TOL = 1e-5


def chrom_counts(prot: pl.DataFrame, vi_unique: pl.DataFrame) -> dict[str, int]:
    m = prot.join(vi_unique, left_on=["Pos", "effectAllele"], right_on=["Pos", "A1"], how="left").with_columns(
        ok=(pl.min_horizontal("ImpFreqA1", 1 - pl.col("ImpFreqA1")) - pl.col("ImpMAF")).abs() <= TOL)
    agreeing_f = pl.col("ImpFreqA1").filter(pl.col("ok"))
    g = m.group_by("r").agg(n=pl.col("ImpFreqA1").is_not_null().sum(), k=pl.col("ok").sum(),
                            f_spread=agreeing_f.max() - agreeing_f.min())
    multi = pl.col("n") > 1
    return g.select(
        rows=pl.len(),
        no_match=(pl.col("n") == 0).sum(),
        single=(pl.col("n") == 1).sum(),
        single_agrees=((pl.col("n") == 1) & (pl.col("k") == 1)).sum(),
        multi=multi.sum(),
        multi_one_agrees=(multi & (pl.col("k") == 1)).sum(),
        multi_none_agree=(multi & (pl.col("k") == 0)).sum(),
        multi_several_agree=(multi & (pl.col("k") > 1)).sum(),
        multi_several_agree_same_f=(multi & (pl.col("k") > 1) & (pl.col("f_spread") <= TOL)).sum(),
    ).row(0, named=True)


def main() -> None:
    cache = Path(sys.argv[1])
    p = pl.scan_parquet(cache / "protein.parquet").with_row_index("r")
    tot: dict[str, int] = {}
    for chrom in p.select(pl.col("Chrom").unique()).collect()["Chrom"]:
        prot = p.filter(pl.col("Chrom") == chrom).drop("Name").collect()
        vi = (pl.scan_parquet(cache / "variant_info.parquet").filter(pl.col("Chrom") == chrom)
              .drop("Chrom").unique().collect())
        for k, x in chrom_counts(prot, vi).items():
            tot[k] = tot.get(k, 0) + x
    for k, x in tot.items():
        print(f"{k:28s} {x:>12,} {100 * x / tot['rows']:8.4f}%")


if __name__ == "__main__":
    main()
