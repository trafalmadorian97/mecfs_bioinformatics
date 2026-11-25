"""
What is in the Pan-UKBB variant manifest, and is it aligned with the hg19 FASTA?

Throwaway probe for the AF reference brainstorm. The manifest is
https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz
(documented as GRCh37, ref = forward-strand reference allele). This script reports row
counts per contig, allele classes, per-population AF availability and AN behaviour, and
compares every ref against the UCSC hg19 FASTA, classifying mismatches (allele swap,
N in FASTA, other).

Usage (on the decompressed manifest):
    pixi r python experiments/claude/pan_ukbb_manifest/inspect_manifest.py MANIFEST.tsv FASTA \
        2>&1 | tee experiments/claude/pan_ukbb_manifest/inspect_manifest.log
"""

import sys
from pathlib import Path

import polars as pl
import pysam

POPS = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID"]


def load(manifest: Path) -> pl.DataFrame:
    columns = (
        ["chrom", "pos", "ref", "alt", "info", "high_quality", "pass_gnomad_genomes"]
        + [f"af_{p}" for p in POPS]
        + [f"an_{p}" for p in POPS]
        + ["gnomad_genomes_af_EUR"]
    )
    return (
        pl.scan_csv(
            manifest,
            separator="\t",
            null_values=["NA"],
            schema_overrides={"chrom": pl.String},
            infer_schema_length=100_000,
        )
        .select(columns)
        .collect(engine="streaming")
    )


def describe(df: pl.DataFrame) -> None:
    print(f"{df.height} rows")
    print("\n== rows per contig")
    with pl.Config(tbl_rows=30):
        print(df.group_by("chrom").len().sort("chrom"))
    print("\n== duplicate chrom/pos/ref/alt:", df.height - df.unique(["chrom", "pos", "ref", "alt"]).height)
    print("\n== allele classes")
    print(
        df.select(
            snv=((pl.col("ref").str.len_bytes() == 1) & (pl.col("alt").str.len_bytes() == 1)).sum(),
            indel=((pl.col("ref").str.len_bytes() != 1) | (pl.col("alt").str.len_bytes() != 1)).sum(),
            non_acgt_ref=(~pl.col("ref").str.contains(r"^[ACGT]+$")).sum(),
            non_acgt_alt=(~pl.col("alt").str.contains(r"^[ACGT]+$")).sum(),
        )
    )
    print("\n== info")
    print(df.select(pl.col("info").min().alias("min"), pl.col("info").median().alias("median"), pl.col("info").null_count().alias("nulls")))
    print("\n== per-population af non-null count, an distinct values, an range")
    for p in POPS:
        print(
            p,
            df[f"af_{p}"].is_not_null().sum(),
            df[f"an_{p}"].n_unique(),
            df[f"an_{p}"].min(),
            df[f"an_{p}"].max(),
            "af==0:",
            (df[f"af_{p}"] == 0).sum(),
        )
    print("\n== high_quality / pass_gnomad_genomes")
    print(df.group_by("high_quality", "pass_gnomad_genomes").len().sort("len", descending=True))
    print("\n== EUR maf distribution")
    maf = pl.min_horizontal(pl.col("af_EUR"), 1 - pl.col("af_EUR"))
    print(
        df.select(
            lt_1e3=(maf < 1e-3).sum(),
            lt_1e2=(maf < 1e-2).sum(),
            ge_1e2=(maf >= 1e-2).sum(),
        )
    )


def fasta_check(df: pl.DataFrame, fasta_path: Path) -> None:
    print("\n== ref vs FASTA")
    results = []
    with pysam.FastaFile(str(fasta_path)) as fasta:
        for (chrom,), part in df.select("chrom", "pos", "ref", "alt").partition_by("chrom", as_dict=True).items():
            contig = f"chr{chrom}"
            seq = fasta.fetch(contig).upper()
            pos = part["pos"].to_list()
            ref = part["ref"].to_list()
            alt = part["alt"].to_list()
            fasta_ref = [seq[p - 1 : p - 1 + len(r)] for p, r in zip(pos, ref)]
            fasta_alt = [seq[p - 1 : p - 1 + len(a)] for p, a in zip(pos, alt)]
            results.append(
                part.with_columns(
                    pl.Series("fasta_ref", fasta_ref),
                    pl.Series("fasta_alt", fasta_alt),
                )
            )
    checked = pl.concat(results).with_columns(
        ref_match=pl.col("fasta_ref") == pl.col("ref"),
        alt_match=pl.col("fasta_alt") == pl.col("alt"),
        fasta_has_n=pl.col("fasta_ref").str.contains("N"),
    )
    print(checked.group_by("ref_match", "alt_match", "fasta_has_n").len().sort("len", descending=True))
    mismatches = checked.filter(~pl.col("ref_match"))
    print(f"{mismatches.height} ref mismatches; by contig:")
    print(mismatches.group_by("chrom").len().sort("chrom"))
    with pl.Config(tbl_rows=30, fmt_str_lengths=30):
        print(mismatches.head(30))


def main(manifest: Path, fasta_path: Path) -> None:
    df = load(manifest)
    describe(df)
    fasta_check(df, fasta_path)


if __name__ == "__main__":
    assert len(sys.argv) == 3, "usage: inspect_manifest.py MANIFEST.tsv FASTA"
    main(Path(sys.argv[1]), Path(sys.argv[2]))
