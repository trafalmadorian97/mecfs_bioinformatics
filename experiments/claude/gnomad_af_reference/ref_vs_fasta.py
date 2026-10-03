"""
Do gnomAD REF alleles match our reference FASTA exactly?

Throwaway probe for the gnomAD AF reference brainstorm, deciding whether the extraction
Task can assert zero REF/FASTA mismatches. It compares every REF in a chr21 panel written
by measure_chr21_extract.py with the FASTA sequence at POS..POS+len(REF)-1 (uppercased,
since UCSC soft-masks repeats in lowercase), and classifies mismatches by what the FASTA
holds there: only ACGT, any N, or another IUPAC ambiguity code. A mismatch over pure ACGT
reference sequence has no benign explanation.

Usage:
    pixi r python experiments/claude/gnomad_af_reference/ref_vs_fasta.py PANEL.parquet FASTA CONTIG \
        2>&1 | tee experiments/claude/gnomad_af_reference/ref_vs_fasta_RELEASE.log
"""

import sys
from pathlib import Path

import polars as pl
import pysam

_ACGT = set("ACGT")


def fasta_class(bases: str) -> str:
    if set(bases) <= _ACGT:
        return "ACGT"
    if "N" in bases:
        return "contains N"
    return "other IUPAC"


def main(panel_path: Path, fasta_path: Path, contig: str) -> None:
    panel = pl.read_parquet(panel_path, columns=["POS", "REF", "ALT"])
    with pysam.FastaFile(str(fasta_path)) as fasta:
        sequence = fasta.fetch(contig).upper()
    print(f"{panel.height} records; {contig} length {len(sequence)}")
    observed = [
        sequence[pos - 1 : pos - 1 + len(ref)]
        for pos, ref in zip(panel["POS"].to_list(), panel["REF"].to_list())
    ]
    checked = panel.with_columns(
        pl.Series("FASTA", observed),
        pl.col("REF").str.len_bytes().alias("ref_len"),
    ).with_columns(
        pl.col("FASTA")
        .map_elements(fasta_class, return_dtype=pl.String)
        .alias("fasta_class"),
        (pl.col("FASTA") != pl.col("REF")).alias("mismatch"),
    )
    print("\n== FASTA class of every REF span")
    print(checked.group_by("fasta_class").agg(pl.len(), pl.col("mismatch").sum()))
    mismatches = checked.filter(pl.col("mismatch"))
    print(f"\n== {mismatches.height} mismatches")
    with pl.Config(tbl_rows=30, fmt_str_lengths=40):
        print(mismatches.head(30))
    print("\n== non-ACGT characters in REF / ALT")
    for column in ["REF", "ALT"]:
        bad = panel.filter(~pl.col(column).str.contains(r"^[ACGT]+$"))
        print(column, bad.height, bad[column].unique().head(10).to_list())


if __name__ == "__main__":
    assert len(sys.argv) == 4, "usage: ref_vs_fasta.py PANEL.parquet FASTA CONTIG"
    main(Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3])
