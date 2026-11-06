"""
Probe: does the non-effect allele of UKB-PPP (regenie ALLELE0) and Western CSF
(GWAS-SSF other_allele) equal the hg38 FASTA reference base at the reported position?

Motivation: the planned large-scale colocalization pipeline joins trait and protein
variants exactly on (chr, pos, EA, NEA) in GRCh38, with NEA == hg38 REF on both sides.
This measures how often each proteomics source already satisfies that.

For each variant the probe classifies (nea_is_ref, ea_is_ref):
  - nea_only: NEA matches REF and EA does not (the desired orientation)
  - ea_only: EA matches REF and NEA does not (swapped; would need a flip)
  - both: both match (expected for indels sharing an anchor base, rare for SNVs)
  - neither: no allele matches (wrong build/position or REF changed)
SNVs (both alleles length 1) and indels are reported separately.

Usage:
  pixi r python experiments/claude/coloc_allele_orientation/probe.py \
      --ppp-dir <extracted protein dir> [--ppp-dir ...] \
      --csf-file <GCST.tsv.gz> [--csf-file ...] 2>&1 | tee experiments/claude/logs/coloc_allele_orientation_probe.log
"""

import argparse
from pathlib import Path

import numpy as np
import polars as pl

from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
    reference_matches,
)

HG38_FASTA_DIR = Path(
    "assets/base_asset_store/reference_data/genome_sequence/ucsc_hg38/processed/"
    "ucsc_hg38_indexed_fasta_wrapped_remove_deps"
)
HG19_FASTA_DIR = Path(
    "assets/base_asset_store/reference_data/genome_sequence/ucsc_hg19/processed/"
    "ucsc_hg19_indexed_fasta_wrapped_remove_deps"
)
CATEGORIES =["nea_only", "ea_only", "both", "neither"]


def classify(fasta: IndexedFasta, frame: pl.DataFrame) -> pl.DataFrame:
    """Per-variant category, given columns chrom (gwaslab code), pos, ea, nea."""
    parts = []
    for (chrom,), sub in frame.group_by(["chrom"]):
        assert isinstance(chrom, int)
        positions = sub["pos"].to_numpy()
        nea_ref = reference_matches(fasta, chrom, positions, sub["nea"])
        ea_ref = reference_matches(fasta, chrom, positions, sub["ea"])
        category = np.where(
            nea_ref & ~ea_ref,
            "nea_only",
            np.where(ea_ref & ~nea_ref, "ea_only", np.where(nea_ref, "both", "neither")),
        )
        parts.append(
            sub.with_columns(category=pl.Series(category)).select(
                "chrom", "pos", "ea", "nea", "category"
            )
        )
    return pl.concat(parts).with_columns(
        is_snv=(pl.col("ea").str.len_bytes() == 1) & (pl.col("nea").str.len_bytes() == 1)
    )


def summarize(label: str, classified: pl.DataFrame) -> None:
    print(f"\n=== {label}: {classified.height:,} variants ===")
    for is_snv, kind in [(True, "SNV"), (False, "indel")]:
        sub = classified.filter(pl.col("is_snv") == is_snv)
        if sub.height == 0:
            continue
        counts = dict(sub["category"].value_counts().iter_rows())
        cells = ", ".join(
            f"{c}={counts.get(c, 0):,} ({100 * counts.get(c, 0) / sub.height:.4f}%)"
            for c in CATEGORIES
        )
        print(f"  {kind} (n={sub.height:,}): {cells}")
    bad = classified.filter(pl.col("is_snv") & pl.col("category").is_in(["ea_only", "neither"]))
    if bad.height:
        print("  example non-conforming SNVs:")
        print(bad.head(8))


def load_ppp(protein_dir: Path) -> pl.DataFrame:
    """All chromosomes of one extracted regenie protein directory."""
    frames = []
    for path in sorted(protein_dir.glob("*.gz")):
        frames.append(
            pl.read_csv(
                path,
                separator=" ",
                columns=["CHROM", "GENPOS", "ALLELE0", "ALLELE1"],
                schema_overrides={"CHROM": pl.Int64, "GENPOS": pl.Int64},
            ).select(
                chrom=pl.col("CHROM"),
                pos=pl.col("GENPOS"),
                ea=pl.col("ALLELE1"),
                nea=pl.col("ALLELE0"),
            )
        )
    assert frames, f"no .gz files in {protein_dir}"
    return pl.concat(frames)


def ppp_hg19_crosscheck(
    hg19: IndexedFasta, protein_dir: Path, classified: pl.DataFrame
) -> None:
    """For hg38-non-conforming SNVs, is ALLELE0 the hg19 REF at the ID's hg19 position?

    If so, ALLELE0 is consistently the hg19 reference allele, and the hg38 mismatch comes
    from the reference base or position changing between builds.
    """
    bad = classified.filter(
        pl.col("is_snv") & pl.col("category").is_in(["ea_only", "neither"])
    ).select("chrom", "pos", "category")
    ids = pl.concat(
        [
            pl.read_csv(
                path,
                separator=" ",
                columns=["CHROM", "GENPOS", "ID"],
                schema_overrides={"CHROM": pl.Int64, "GENPOS": pl.Int64},
            )
            for path in sorted(protein_dir.glob("*.gz"))
        ]
    ).rename({"CHROM": "chrom", "GENPOS": "pos"})
    parsed = (
        bad.join(ids, on=["chrom", "pos"], how="left")
        .unique(["chrom", "pos", "ID"])
        .with_columns(parts=pl.col("ID").str.split(":"))
        .with_columns(
            pos19=pl.col("parts").list.get(1).cast(pl.Int64),
            a0=pl.col("parts").list.get(2),
            a1=pl.col("parts").list.get(3),
        )
        .filter((pl.col("a0").str.len_bytes() == 1) & (pl.col("a1").str.len_bytes() == 1))
    )
    rows = classify(
        hg19, parsed.select("chrom", pos="pos19", ea="a1", nea="a0")
    ).rename({"category": "hg19_category"})
    print(f"  hg19 crosscheck of {rows.height:,} hg38-non-conforming SNVs (ID position):")
    joined = parsed.select("chrom", "pos19", "category").join(
        rows.select("chrom", pos19="pos", hg19_category="hg19_category"),
        on=["chrom", "pos19"],
        how="left",
    )
    print(
        joined.group_by("category", "hg19_category")
        .len()
        .sort("category", "hg19_category")
    )


def load_csf(path: Path) -> pl.DataFrame:
    """One GWAS-SSF file (autosomes only)."""
    return pl.read_csv(
        path,
        separator="\t",
        columns=["chromosome", "base_pair_location", "effect_allele", "other_allele"],
        schema_overrides={"chromosome": pl.Int64, "base_pair_location": pl.Int64},
    ).select(
        chrom=pl.col("chromosome"),
        pos=pl.col("base_pair_location"),
        ea=pl.col("effect_allele"),
        nea=pl.col("other_allele"),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ppp-dir", type=Path, action="append", default=[])
    parser.add_argument("--csf-file", type=Path, action="append", default=[])
    args = parser.parse_args()
    fasta = IndexedFasta.open(HG38_FASTA_DIR)
    hg19 = IndexedFasta.open(HG19_FASTA_DIR)
    for protein_dir in args.ppp_dir:
        classified = classify(fasta, load_ppp(protein_dir))
        summarize(f"PPP {protein_dir.name}", classified)
        ppp_hg19_crosscheck(hg19, protein_dir, classified)
    for path in args.csf_file:
        summarize(f"CSF {path.name}", classify(fasta, load_csf(path)))


if __name__ == "__main__":
    main()
