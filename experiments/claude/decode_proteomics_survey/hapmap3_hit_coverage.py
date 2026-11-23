"""
How many genome-wide-significant protein loci would a HapMap3-only screen miss?

The deCODE hit table (like PPP's and CSF's) finds protein hit segments from the shared HapMap3
index only: significant index variants (p < 5e-8) merged when closer than 250 kb. A locus
whose signal no HapMap3 variant picks up (e.g. a rare, Iceland-enriched variant with few
common proxies) never becomes a coloc candidate.

For each protein file this script:
  1. reads the autosomal rows with p < 5e-8 (by minus_log10_pval);
  2. clusters them into loci the way the hit table does (merge gap 250 kb);
  3. marks which significant rows are HapMap3 variants (gwaslab's hg38 snplist, joined on
     chromosome, position and the unordered allele pair, as the index is built);
  4. reports loci with no significant HapMap3 variant ("missed"), with the lead's p and ImpMAF.

Each protein's loci are appended to <loci.csv> (one row per locus) for later summarizing
with summarize_hapmap3_hit_coverage.py.

Usage (from the repo root):
  pixi r python experiments/claude/decode_proteomics_survey/hapmap3_hit_coverage.py \
      <loci.csv> <protein.txt.gz> [...] | tee -a experiments/claude/decode_proteomics_survey/hapmap3_hit_coverage.log
"""

import math
import sys
from pathlib import Path

import polars as pl

from mecfs_bio.build_system.task.ppp_database.hapmap3_membership_task import (
    bundled_hapmap3_hg38_snplist_path,
    normalize_hapmap3_snplist,
)

MIN_MLOG10P = -math.log10(5e-8)
MERGE_GAP_BP = 250_000


def allele_pair(a: str, b: str) -> list[pl.Expr]:
    return [
        pl.min_horizontal(pl.col(a), pl.col(b)).alias("_a1"),
        pl.max_horizontal(pl.col(a), pl.col(b)).alias("_a2"),
    ]


def hapmap3_keys() -> pl.DataFrame:
    m = normalize_hapmap3_snplist(bundled_hapmap3_hg38_snplist_path())
    m = m if isinstance(m, pl.DataFrame) else m.collect()
    return (
        m.select(pl.col("CHR").cast(pl.Int64), pl.col("POS").cast(pl.Int64), *allele_pair("EA", "NEA"))
        .unique()
        .with_columns(in_hm3=pl.lit(True))
    )


def significant_rows(path: Path, hm3: pl.DataFrame) -> pl.DataFrame:
    rows = (
        pl.scan_csv(path, separator="\t", schema_overrides={"minus_log10_pval": pl.Float64, "Pos": pl.Int64})
        .filter((pl.col("minus_log10_pval") >= MIN_MLOG10P) & (pl.col("Chrom") != "chrX"))
        .select(
            pl.col("Chrom").str.strip_prefix("chr").cast(pl.Int64).alias("CHR"),
            pl.col("Pos").alias("POS"),
            "effectAllele", "otherAllele", "ImpMAF", "minus_log10_pval",
        )
        .collect()
        .with_columns(allele_pair("effectAllele", "otherAllele"))
    )
    return (
        rows.join(hm3, on=["CHR", "POS", "_a1", "_a2"], how="left")
        .with_columns(pl.col("in_hm3").fill_null(False))
        .sort("CHR", "POS")
    )


def loci(rows: pl.DataFrame) -> pl.DataFrame:
    new_locus = (pl.col("CHR") != pl.col("CHR").shift()) | (
        pl.col("POS") - pl.col("POS").shift() > MERGE_GAP_BP
    )
    with_id = rows.with_columns(locus=new_locus.fill_null(True).cum_sum())
    lead = with_id.sort("minus_log10_pval", descending=True).group_by("locus", maintain_order=True).first()
    summary = with_id.group_by("locus").agg(
        start=pl.col("POS").min(),
        end=pl.col("POS").max(),
        n_sig=pl.len(),
        n_sig_hm3=pl.col("in_hm3").sum(),
        best_hm3_mlog10p=pl.col("minus_log10_pval").filter(pl.col("in_hm3")).max(),
    )
    return summary.join(
        lead.select("locus", "CHR", lead_pos="POS", lead_mlog10p="minus_log10_pval", lead_maf="ImpMAF",
                    lead_in_hm3="in_hm3"),
        on="locus",
    ).sort("CHR", "start")


def main() -> None:
    hm3 = hapmap3_keys()
    pl.Config.set_tbl_rows(50)
    pl.Config.set_tbl_cols(20)
    totals = {"proteins": 0, "loci": 0, "missed": 0}
    out_csv = Path(sys.argv[1])
    for path in (Path(a) for a in sys.argv[2:]):
        found = loci(significant_rows(path, hm3))
        rows = found.drop("locus").with_columns(protein=pl.lit(path.name))
        with open(out_csv, "a") as out:
            rows.write_csv(out, include_header=out.tell() == 0)
        missed = found.filter(pl.col("n_sig_hm3") == 0)
        totals["proteins"] += 1
        totals["loci"] += found.height
        totals["missed"] += missed.height
        print(f"\n{path.name}: {found.height} autosomal loci at p<5e-8, {missed.height} with no "
              f"significant HapMap3 variant")
        print(found.drop("locus"))
    print(f"\nTOTAL {totals}")


if __name__ == "__main__":
    main()
