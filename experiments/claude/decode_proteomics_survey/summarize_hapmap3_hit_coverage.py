"""
Summarize the per-locus CSV written by hapmap3_hit_coverage.py: how many protein loci at
p < 5e-8 have no significant HapMap3 variant, broken down by the lead variant's frequency,
by the number of significant variants in the locus, and by the lead's strength.

Usage (from the repo root):
  pixi r python experiments/claude/decode_proteomics_survey/summarize_hapmap3_hit_coverage.py <loci.csv>
"""

import sys

import polars as pl


def table(loci: pl.DataFrame, by: pl.Expr) -> pl.DataFrame:
    return (
        loci.group_by(by.alias("bin"))
        .agg(loci=pl.len(), missed=(pl.col("n_sig_hm3") == 0).sum())
        .with_columns(missed_pct=(100 * pl.col("missed") / pl.col("loci")).round(1))
        .sort("bin")
    )


def main() -> None:
    loci = pl.read_csv(sys.argv[1])
    pl.Config.set_tbl_rows(100)
    print(f"proteins: {loci['protein'].n_unique()}, loci: {loci.height}, "
          f"missed: {(loci['n_sig_hm3'] == 0).sum()}")
    maf = pl.col("lead_maf")
    print("\nby lead MAF")
    print(table(loci, pl.when(maf < 0.001).then(pl.lit("1: <0.1%"))
                .when(maf < 0.01).then(pl.lit("2: 0.1-1%"))
                .when(maf < 0.05).then(pl.lit("3: 1-5%"))
                .otherwise(pl.lit("4: >=5%"))))
    n = pl.col("n_sig")
    print("\nby number of significant variants in the locus")
    print(table(loci, pl.when(n == 1).then(pl.lit("1: 1"))
                .when(n <= 5).then(pl.lit("2: 2-5"))
                .when(n <= 50).then(pl.lit("3: 6-50"))
                .otherwise(pl.lit("4: >50"))))
    p = pl.col("lead_mlog10p")
    print("\nby lead -log10 p")
    print(table(loci, pl.when(p < 10).then(pl.lit("1: 7.3-10"))
                .when(p < 20).then(pl.lit("2: 10-20"))
                .when(p < 50).then(pl.lit("3: 20-50"))
                .otherwise(pl.lit("4: >=50"))))
    print("\nmissed loci with >= 6 significant variants or lead -log10 p >= 20 (more credible)")
    credible = loci.filter((pl.col("n_sig_hm3") == 0) & ((n >= 6) | (p >= 20)))
    print(credible.select("protein", "CHR", "lead_pos", "n_sig", "lead_mlog10p", "lead_maf")
          .with_columns(pl.col("protein").str.extract(r"PC0_\d+_\d+_([^_]+)_")).sort("lead_mlog10p", descending=True))
    print("\nlargest -log10 p of any missed lead per protein vs its strongest locus")
    print(loci.group_by("protein").agg(
        strongest=pl.col("lead_mlog10p").max(),
        strongest_missed=pl.col("lead_mlog10p").filter(pl.col("n_sig_hm3") == 0).max(),
    ).with_columns(pl.col("protein").str.extract(r"PC0_\d+_\d+_([^_]+)_")).sort("strongest_missed", descending=True))


if __name__ == "__main__":
    main()
