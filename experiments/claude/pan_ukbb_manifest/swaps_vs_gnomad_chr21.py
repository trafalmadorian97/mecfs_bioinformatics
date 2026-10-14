"""
Which allele does a swapped Pan-UKBB row's af describe? Check chr21 against gnomAD v2.1.1.

Takes the Pan-UKBB chr21 rows whose ref disagrees with the hg19 FASTA and looks up every
gnomAD v2.1.1 genome record (all FILTERs) at the same position, from the chr21 TSV written
by experiments/claude/gnomad_af_reference/measure_chr21_extract.py. For a pure label swap
(the manifest's alt is really the reference base, and af describes it), af_EUR should be
close to 1 - gnomAD AF_nfe of the record whose ALT is the manifest's ref. If instead af
describes the non-reference allele, af_EUR should be close to gnomAD AF_nfe itself.

Usage:
    pixi r python experiments/claude/pan_ukbb_manifest/swaps_vs_gnomad_chr21.py MANIFEST.tsv FASTA GNOMAD_CHR21_TSV \
        2>&1 | tee experiments/claude/pan_ukbb_manifest/swaps_vs_gnomad_chr21.log
"""

import sys
from pathlib import Path

import polars as pl
import pysam

GNOMAD_GROUPS = ["afr", "amr", "asj", "eas", "fin", "nfe", "oth", "nfe_nwe", "nfe_seu", "nfe_onf", "nfe_est"]


def swapped_chr21(manifest: Path, fasta_path: Path) -> pl.DataFrame:
    rows = (
        pl.scan_csv(manifest, separator="\t", null_values=["NA"], schema_overrides={"chrom": pl.String}, infer_schema_length=100_000)
        .filter(pl.col("chrom") == "21")
        .select("pos", "ref", "alt", "rsid", "af_EUR")
        .collect()
    )
    with pysam.FastaFile(str(fasta_path)) as fasta:
        seq = fasta.fetch("chr21").upper()
    fasta_ref = [seq[p - 1 : p - 1 + len(r)] for p, r in zip(rows["pos"].to_list(), rows["ref"].to_list())]
    return rows.with_columns(pl.Series("fasta_ref", fasta_ref)).filter(pl.col("fasta_ref") != pl.col("ref"))


def gnomad_chr21(path: Path) -> pl.DataFrame:
    names = ["CHROM", "POS", "REF", "ALT", "FILTER"]
    for g in GNOMAD_GROUPS:
        names += [f"AF_{g}", f"AN_{g}"]
    return pl.read_csv(
        path, separator="\t", has_header=False, new_columns=names, null_values=["."],
        schema_overrides={"CHROM": pl.String, "POS": pl.Int64, "FILTER": pl.String},
        infer_schema_length=0,
    ).select("POS", "REF", "ALT", "FILTER", pl.col("AF_nfe").cast(pl.Float64), pl.col("AN_nfe").cast(pl.Int64))


def main(manifest: Path, fasta_path: Path, gnomad_tsv: Path) -> None:
    swapped = swapped_chr21(manifest, fasta_path)
    gnomad = gnomad_chr21(gnomad_tsv)
    joined = swapped.join(gnomad, left_on="pos", right_on="POS", how="left")
    print(f"{swapped.height} swapped chr21 rows; {joined['REF'].is_not_null().sum()} gnomAD records at those positions")
    with pl.Config(tbl_rows=100, tbl_cols=20, fmt_str_lengths=12, tbl_width_chars=200):
        print(joined.sort("pos"))
    label_swap = joined.filter((pl.col("REF") == pl.col("alt")) & (pl.col("ALT") == pl.col("ref")))
    print(f"\n== {label_swap.height} rows where gnomAD has the same alleles in FASTA orientation")
    print(
        label_swap.select(
            n=pl.len(),
            af_eur_close_to_1_minus_gnomad=((pl.col("af_EUR") - (1 - pl.col("AF_nfe"))).abs() < 0.05).sum(),
            af_eur_close_to_gnomad=((pl.col("af_EUR") - pl.col("AF_nfe")).abs() < 0.05).sum(),
            gnomad_pass=(pl.col("FILTER") == "PASS").sum(),
        )
    )


if __name__ == "__main__":
    assert len(sys.argv) == 4, "usage: swaps_vs_gnomad_chr21.py MANIFEST.tsv FASTA GNOMAD_CHR21_TSV"
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
