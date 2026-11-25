"""
(1) Structure of variant_info_iceland per position: how many entries, whether the GRCh38
    reference allele is among them, and whether the entries sum to 1 (within 1e-3).
(2) The strict cleanliness rule for deCODE rows (reference file AP4M1): a row is a clean
    two-allele test when variant_info lists nothing at its position except the row's own
    effect allele. Share of rows that pass, per row kind, alone vs at multi-row positions, MAF.

  pixi r python experiments/claude/decode_proteomics_survey/variant_info_record_structure.py
"""

import numpy as np
import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.decode_somascan.raw.decode_reference_protein_file import (
    DECODE_SMP_AP4M1_RAW_FILE,
)
from mecfs_bio.assets.gwas.decode_somascan.raw.decode_variant_info import (
    DECODE_VARIANT_INFO,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.decode_database.decode_sumstats_format import (
    chrom_code_expr,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
    reference_matches,
)

pl.Config.set_tbl_rows(40)
pl.Config.set_tbl_width_chars(220)
SUM_TOLERANCE = 1e-3


def with_is_ref(frame: pl.DataFrame, fasta: IndexedFasta, allele_col: str) -> pl.DataFrame:
    parts = []
    for (chrom,), part in frame.group_by("CHR", maintain_order=True):
        alleles = part[allele_col]
        ok = (alleles != "*").to_numpy()
        is_ref = np.zeros(part.height, dtype=bool)
        is_ref[ok] = reference_matches(
            fasta, int(chrom), part["POS"].to_numpy()[ok], alleles.filter(ok)
        )
        parts.append(part.with_columns(pl.Series("is_ref", is_ref)))
    return pl.concat(parts)


def structure(vi: pl.DataFrame) -> pl.DataFrame:
    per_pos = vi.group_by("CHR", "POS").agg(
        pl.len().alias("n_entries"),
        pl.col("is_ref").any().alias("ref_listed"),
        ((pl.col("ImpFreqA1").sum() - 1).abs() <= SUM_TOLERANCE).alias("sums_to_1"),
        (pl.col("A1") == "*").any().alias("star_listed"),
    )
    print("\n(1) variant_info positions by entry count and reference listing")
    print(
        per_pos.with_columns(
            pl.when(pl.col("n_entries") >= 4).then(pl.lit("4+")).otherwise(
                pl.col("n_entries").cast(pl.String)
            ).alias("entries")
        )
        .group_by("entries", "ref_listed")
        .agg(
            pl.len().alias("positions"),
            pl.col("sums_to_1").mean().alias("share_sum_to_1"),
            pl.col("star_listed").mean().alias("share_with_star"),
        )
        .sort("entries", "ref_listed")
    )
    return per_pos


def strict_rule(rows: pl.DataFrame, vi: pl.DataFrame) -> None:
    only = (
        vi.group_by("CHR", "POS")
        .agg(pl.len().alias("n_entries"), pl.col("A1").first().alias("only_allele"))
    )
    rows = rows.join(only, on=["CHR", "POS"], how="left").with_columns(
        ((pl.col("n_entries") == 1) & (pl.col("only_allele") == pl.col("EA")))
        .fill_null(False)
        .alias("clean")
    )
    print("\n(2) strict rule: variant_info lists only the row's own allele at its position")
    print(
        rows.group_by("kind", "lone", "common")
        .agg(pl.len().alias("rows"), pl.col("clean").sum().alias("clean"), pl.col("clean").mean().alias("share_clean"))
        .sort("kind", "lone", "common", descending=[False, True, True])
    )
    for kind in ("snv", "indel"):
        k = rows.filter(pl.col("kind") == kind, "common")
        print(
            f"common {kind}s kept by the strict rule: {int(k['clean'].sum())} of {k.height} "
            f"({k['clean'].mean():.1%}); kept by the multi-row rule: "
            f"{int(k['lone'].sum())} ({k['lone'].mean():.1%})"
        )


def main() -> None:
    assets = DEFAULT_RUNNER.run(
        [DECODE_SMP_AP4M1_RAW_FILE, DECODE_VARIANT_INFO, UCSC_HG38_INDEXED_FASTA]
    )
    raw = assets[DECODE_SMP_AP4M1_RAW_FILE.asset_id]
    vi_asset = assets[DECODE_VARIANT_INFO.asset_id]
    fasta_dir = assets[UCSC_HG38_INDEXED_FASTA.asset_id]
    assert isinstance(raw, FileAsset) and isinstance(vi_asset, FileAsset)
    assert isinstance(fasta_dir, DirectoryAsset)
    fasta = IndexedFasta.open(fasta_dir.path)

    vi = with_is_ref(
        pl.read_parquet(vi_asset.path, columns=["CHR", "POS", "A1", "ImpFreqA1"]), fasta, "A1"
    )
    structure(vi)

    ea, oa = pl.col("effectAllele"), pl.col("otherAllele")
    rows = (
        pl.read_csv(
            raw.path,
            separator="\t",
            columns=["Chrom", "Pos", "effectAllele", "otherAllele", "ImpMAF"],
            null_values=["."],
            quote_char=None,
            schema_overrides={"Pos": pl.Int64},
        )
        .with_columns(pl.len().over("Chrom", "Pos").alias("n_at_pos"))
        .filter((ea != "*") & (oa != "*") & (ea != oa))
        .select(
            chrom_code_expr("Chrom").alias("CHR"),
            pl.col("Pos").alias("POS"),
            ea.alias("EA"),
            pl.when((ea.str.len_bytes() == 1) & (oa.str.len_bytes() == 1))
            .then(pl.lit("snv"))
            .otherwise(pl.lit("indel"))
            .alias("kind"),
            (pl.col("n_at_pos") == 1).alias("lone"),
            (pl.col("ImpMAF") >= 0.01).fill_null(False).alias("common"),
        )
    )
    strict_rule(rows, vi)


if __name__ == "__main__":
    main()
