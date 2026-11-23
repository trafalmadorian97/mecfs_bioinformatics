"""
Count the rows of one deCODE SomaScan sumstats file at or below a few p-value thresholds, on
chrX and over all chromosomes. This sizes the per-protein chrX supplement (p <= 1e-5) in the
database design.

Filters on minus_log10_pval, not Pval: Pval rounds to 0 for very strong hits.

Usage (from the repo root), on a downloaded file:
  pixi r python experiments/claude/decode_proteomics_survey/count_significant_rows.py \
      <protein.txt.gz> | tee experiments/claude/decode_proteomics_survey/count_significant_rows.log
"""

import math
import sys
import time
from pathlib import Path

import polars as pl

THRESHOLDS_MLOG10P = {"p<=1e-5": 5.0, "p<=1e-6": 6.0, "p<=5e-8": -math.log10(5e-8)}


def main() -> None:
    path = Path(sys.argv[1])
    t0 = time.perf_counter()
    frame = (
        pl.scan_csv(path, separator="\t", schema_overrides={"minus_log10_pval": pl.Float64})
        .select(
            is_x=pl.col("Chrom") == "chrX",
            mlog10p=pl.col("minus_log10_pval"),
        )
        .collect()
    )
    print(f"file: {path.name}")
    print(f"read in {time.perf_counter() - t0:.1f}s")
    print(f"{'':10s} {'all chroms':>12s} {'chrX':>10s}")
    print(f"{'rows':10s} {frame.height:>12,} {frame['is_x'].sum():>10,}")
    for label, cut in THRESHOLDS_MLOG10P.items():
        hits = frame.filter(pl.col("mlog10p") >= cut)
        print(f"{label:10s} {hits.height:>12,} {hits['is_x'].sum():>10,}")
    expected = frame.filter(pl.col("is_x")).height * 1e-5
    print(f"chrX rows expected at p<=1e-5 with no association at all: ~{expected:.0f}")


if __name__ == "__main__":
    main()
