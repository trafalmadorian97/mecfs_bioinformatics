"""
Follow-up to inspect_manifest.py: characterise the ref/alt-swapped manifest rows.

Reports AN per contig, and for rows whose ref disagrees with the hg19 FASTA (but alt
matches it) shows the EUR frequency next to gnomAD's, whether the variant also appears in
the manifest in the FASTA orientation, and the rsid, so we can tell whether these are
label swaps, frequency swaps, or a second record of the same site.

Usage:
    pixi r python experiments/claude/pan_ukbb_manifest/inspect_swaps.py MANIFEST.tsv FASTA \
        2>&1 | tee experiments/claude/pan_ukbb_manifest/inspect_swaps.log
"""

import sys
from pathlib import Path

import polars as pl
import pysam

SWAP_CONTIGS = ["21", "22", "X"]


def main(manifest: Path, fasta_path: Path) -> None:
    lazy = pl.scan_csv(
        manifest,
        separator="\t",
        null_values=["NA"],
        schema_overrides={"chrom": pl.String},
        infer_schema_length=100_000,
    )
    print("== AN_EUR / AN_AFR per contig (distinct values)")
    with pl.Config(tbl_rows=30):
        print(
            lazy.group_by("chrom")
            .agg(pl.col("an_EUR").unique().sort(), pl.col("an_AFR").unique().sort())
            .sort("chrom")
            .collect(engine="streaming")
        )
    df = (
        lazy.filter(pl.col("chrom").is_in(SWAP_CONTIGS))
        .select(
            "chrom", "pos", "ref", "alt", "rsid", "info", "af_EUR", "gnomad_genomes_af_EUR",
            "pass_gnomad_genomes", "high_quality",
        )
        .collect(engine="streaming")
    )
    parts = []
    with pysam.FastaFile(str(fasta_path)) as fasta:
        for (chrom,), part in df.partition_by("chrom", as_dict=True).items():
            seq = fasta.fetch(f"chr{chrom}").upper()
            parts.append(
                part.with_columns(
                    pl.Series(
                        "fasta_ref",
                        [seq[p - 1 : p - 1 + len(r)] for p, r in zip(part["pos"].to_list(), part["ref"].to_list())],
                    )
                )
            )
    checked = pl.concat(parts)
    swapped = checked.filter(pl.col("fasta_ref") != pl.col("ref"))
    mirror = swapped.join(
        checked.select("chrom", "pos", pl.col("ref").alias("alt"), pl.col("alt").alias("ref"), pl.lit(True).alias("mirror_present")),
        on=["chrom", "pos", "ref", "alt"],
        how="left",
    )
    print(f"\n== {swapped.height} swapped rows")
    print("mirror (FASTA-oriented) record also present:", mirror["mirror_present"].fill_null(False).sum())
    print(swapped.group_by("chrom", "high_quality", "pass_gnomad_genomes").len().sort("chrom"))
    print("gnomad af non-null:", swapped["gnomad_genomes_af_EUR"].is_not_null().sum())
    print("\n== EUR af vs gnomAD af (swapped rows with gnomAD af)")
    comparable = swapped.filter(pl.col("gnomad_genomes_af_EUR").is_not_null())
    print(
        comparable.select(
            n=pl.len(),
            close_same=((pl.col("af_EUR") - pl.col("gnomad_genomes_af_EUR")).abs() < 0.05).sum(),
            close_flipped=((pl.col("af_EUR") - (1 - pl.col("gnomad_genomes_af_EUR"))).abs() < 0.05).sum(),
        )
    )
    print("\n== position distribution of swaps (Mb)")
    print(swapped.group_by("chrom", (pl.col("pos") // 1_000_000).alias("mb")).len().sort("len", descending=True).head(15))
    with pl.Config(tbl_rows=40, tbl_cols=12, fmt_str_lengths=20):
        print(swapped.sample(40, seed=0).sort("chrom", "pos"))


if __name__ == "__main__":
    assert len(sys.argv) == 3, "usage: inspect_swaps.py MANIFEST.tsv FASTA"
    main(Path(sys.argv[1]), Path(sys.argv[2]))
