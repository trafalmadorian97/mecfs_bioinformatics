"""What does an SNV reference-match "trust signal" report for DecodeME before gwaslab harmonization?

Trust signal (proposed for CleanIndelOrientation): over non-palindromic SNVs, the fraction whose NEA
equals the reference base on either strand (NEA == REF, or complement(NEA) == REF). Each SNV is
classified as
  nea_ref          NEA == REF                    (orientation right, strand right)
  nea_ref_revcomp  complement(NEA) == REF        (orientation right, strand wrong)
  ea_ref           EA == REF                     (orientation swapped)
  ea_ref_revcomp   complement(EA) == REF         (orientation swapped, strand wrong)
  neither          no allele matches on either strand

Tables checked (both produced by GWASLabCreateSumstatsTask with basic_check=True, no harmonization):
  - build 38 (minimal processing) against hg38: the source orientation before liftover;
  - build 37 (liftover to 37) against hg19: the input to gwaslab harmonization in the annovar
    rsid-assignment chain.
Variants are matched between the two tables by row order key SNPID to attribute build-37 mismatches
to liftover.

Run:
  pixi r python experiments/claude/polyfun_annot_dedup/decodeme_pre_harmonization_trust_signal.py \
    2>&1 | tee experiments/claude/polyfun_annot_dedup/decodeme_pre_harmonization_trust_signal.log
"""

import gc
import pickle
from pathlib import Path

import numpy as np
import polars as pl

from experiments.claude.polyfun_annot_dedup.benchmark_fasta_ref_check import read_fai, ref_match

D = Path("assets/base_asset_store/gwas/ME_CFS/DecodeME/gwaslab_sumstats")
HG19 = Path.home() / ".gwaslab/hg19.fa"
HG38 = Path.home() / ".gwaslab/hg38.fa"
CHROM_NAMES = {"23": "X", "24": "Y", "25": "M"}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def load(name: str) -> pl.DataFrame:
    with open(D / name, "rb") as f:
        sumstats = pickle.load(f)
    print(f"\n{name}: build={sumstats.meta.get('gwaslab', {}).get('genome_build')} columns={list(sumstats.data.columns)}")
    cols = [c for c in ["SNPID", "CHR", "POS", "EA", "NEA", "EAF", "STATUS"] if c in sumstats.data.columns]
    pdf = sumstats.data[cols].copy()
    del sumstats
    gc.collect()
    for c in ["EA", "NEA", "SNPID"]:
        pdf[c] = pdf[c].astype(str)
    df = pl.from_pandas(pdf).with_columns(pl.col("CHR").cast(pl.Utf8).replace(CHROM_NAMES), pl.col("POS").cast(pl.Int64))
    print(df.head(3))
    return df


def classify(df: pl.DataFrame, fasta: Path) -> pl.DataFrame:
    fai = read_fai(fasta)
    df = df.filter(pl.col("CHR").is_in(list(fai.keys())) & pl.col("POS").is_not_null())
    snv = (pl.col("EA").str.len_bytes() == 1) & (pl.col("NEA").str.len_bytes() == 1)
    acgt = pl.col("EA").is_in(list(COMP)) & pl.col("NEA").is_in(list(COMP))
    df = df.filter(snv & acgt).with_columns(
        pl.col("NEA").replace(COMP).alias("NEA_c"),
        pl.col("EA").replace(COMP).alias("EA_c"),
    )
    df = df.filter(pl.col("EA") != pl.col("NEA_c"))  # drop palindromic
    mm = np.memmap(fasta, dtype=np.uint8, mode="r")
    chrom = df["CHR"].to_numpy()
    pos = df["POS"].to_numpy()
    m = {c: ref_match(mm, fai, chrom, pos, df[c]) for c in ["NEA", "NEA_c", "EA", "EA_c"]}
    cls = np.select(
        [m["NEA"], m["NEA_c"], m["EA"], m["EA_c"]],
        ["nea_ref", "nea_ref_revcomp", "ea_ref", "ea_ref_revcomp"],
        default="neither",
    )
    return df.with_columns(pl.Series("cls", cls))


def report(label: str, df: pl.DataFrame) -> None:
    n = df.height
    counts = df.group_by("cls").len().sort("len", descending=True)
    print(f"\n=== {label}: {n:,} non-palindromic SNVs ===")
    print(counts.with_columns((pl.col("len") / n * 100).round(4).alias("pct")))
    trusted = df.filter(pl.col("cls").is_in(["nea_ref", "nea_ref_revcomp"])).height
    classified = df.filter(pl.col("cls") != "neither").height
    print(f"trust signal (NEA is REF on either strand / classifiable) = {trusted / classified * 100:.4f}%")
    print(f"strict (NEA == REF, same strand) / all = {df.filter(pl.col('cls') == 'nea_ref').height / n * 100:.4f}%")
    by_chr = (
        df.group_by("CHR")
        .agg(pl.len().alias("n"), *[(pl.col("cls") == c).mean().mul(100).round(3).alias(c) for c in ["nea_ref", "nea_ref_revcomp", "ea_ref", "ea_ref_revcomp", "neither"]])
        .sort(pl.col("CHR").cast(pl.Int32, strict=False), nulls_last=True)
    )
    with pl.Config(tbl_rows=30, tbl_cols=-1, tbl_width_chars=200):
        print(by_chr)


b38 = classify(load("decode_me_gwas_1_sumstats_minimal_processing.pickle"), HG38)
report("build 38 vs hg38 (before liftover)", b38)
b38_small = b38.select("SNPID", pl.col("cls").alias("cls38"), pl.col("EA").alias("EA38"), pl.col("NEA").alias("NEA38"))
del b38
gc.collect()

b37 = classify(load("decode_me_gwas_1_sumstats_liftover_to_37.pickle"), HG19)
report("build 37 vs hg19 (after liftover, before gwaslab harmonization)", b37)

j = b37.join(b38_small, on="SNPID", how="inner")
print(f"\n=== build 37 classes cross-tabulated against build 38 classes (joined on SNPID: {j.height:,}) ===")
print(f"alleles unchanged by liftover (EA, NEA identical): {(j['EA'] == j['EA38']).mean() * 100:.4f}% EA, {(j['NEA'] == j['NEA38']).mean() * 100:.4f}% NEA")
with pl.Config(tbl_rows=30):
    print(j.group_by("cls38", "cls").len().sort("len", descending=True))
with pl.Config(tbl_cols=-1, tbl_width_chars=200):
    print(j.filter((pl.col("cls38") == "nea_ref") & (pl.col("cls") != "nea_ref")).head(10))
