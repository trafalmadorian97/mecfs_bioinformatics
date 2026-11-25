"""
Are deCODE SNV rows that are alone at their position truly biallelic? Reference file AP4M1,
genome-wide.

(a) Allele coding: how often is otherAllele (resp. effectAllele) the GRCh38 reference allele,
    per row kind.
(b) Spanning deletions: lone SNVs inside the span of a deletion row of the file (a row whose
    longer allele is the reference sequence: it deletes the bases after the anchor).
(c) variant_info_iceland (every allele genotyped at a position, with its frequency): lone SNVs
    whose position lists exactly their two alleles, each once, with frequencies summing to 1
    within 1e-3; versus extra alleles or a frequency shortfall (pooled "rest").

  pixi r python experiments/claude/decode_proteomics_survey/check_lone_snv_biallelic.py
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

FREQ_TOLERANCE = 1e-3
pl.Config.set_tbl_rows(30)
pl.Config.set_tbl_width_chars(220)


def read_rows(path) -> pl.DataFrame:  # noqa: ANN001
    ea, oa = pl.col("EA"), pl.col("OA")
    return (
        pl.read_csv(
            path,
            separator="\t",
            columns=["Chrom", "Pos", "effectAllele", "otherAllele", "ImpMAF"],
            null_values=["."],
            quote_char=None,
            schema_overrides={"Pos": pl.Int64},
        )
        .rename({"effectAllele": "EA", "otherAllele": "OA"})
        .with_columns(
            chrom_code_expr("Chrom").alias("CHR"),
            pl.when((ea == "*") | (oa == "*"))
            .then(pl.lit("star"))
            .when(ea == oa)
            .then(pl.lit("same_allele"))
            .when((ea.str.len_bytes() == 1) & (oa.str.len_bytes() == 1))
            .then(pl.lit("snv"))
            .otherwise(pl.lit("indel"))
            .alias("kind"),
        )
        .drop("Chrom")
        .with_columns(pl.len().over("CHR", "Pos").alias("n_at_pos"))
    )


def add_reference_matches(rows: pl.DataFrame, fasta: IndexedFasta) -> pl.DataFrame:
    parts = []
    for (chrom,), part in rows.group_by("CHR", maintain_order=True):
        pos = part["Pos"].to_numpy()
        matches = {}
        for col in ("EA", "OA"):
            alleles = part[col]
            ok = (alleles != "*").to_numpy()
            m = np.zeros(part.height, dtype=bool)
            m[ok] = reference_matches(fasta, int(chrom), pos[ok], alleles.filter(ok))
            matches[f"{col}_is_ref"] = m
        parts.append(part.with_columns(**{k: pl.Series(v) for k, v in matches.items()}))
    return pl.concat(parts)


def allele_coding(rows: pl.DataFrame) -> None:
    print("\n(a) allele coding vs GRCh38, per row kind")
    print(
        rows.group_by("kind")
        .agg(
            pl.len().alias("rows"),
            pl.col("OA_is_ref").mean().alias("OA_is_ref"),
            pl.col("EA_is_ref").mean().alias("EA_is_ref"),
            (~pl.col("OA_is_ref") & ~pl.col("EA_is_ref")).mean().alias("neither"),
        )
        .sort("rows", descending=True)
    )


def deletion_cover(rows: pl.DataFrame) -> pl.DataFrame:
    """Per row: True when a deletion row starting earlier on its chromosome spans its position."""
    longer_is_ea = pl.col("EA").str.len_bytes() > pl.col("OA").str.len_bytes()
    deletions = (
        rows.filter(
            (pl.col("kind") == "indel")
            & pl.when(longer_is_ea).then(pl.col("EA_is_ref")).otherwise(pl.col("OA_is_ref"))
        )
        .select(
            "CHR",
            pl.col("Pos").alias("del_pos"),
            (
                pl.col("Pos")
                + pl.max_horizontal(pl.col("EA").str.len_bytes(), pl.col("OA").str.len_bytes())
                - 1
            ).alias("del_end"),
        )
        .sort("CHR", "del_pos")
        .with_columns(pl.col("del_end").cum_max().over("CHR").alias("cover_end"))
    )
    return (
        rows.with_columns((pl.col("Pos") - 1).alias("_before"))
        .sort("CHR", "_before")
        .join_asof(
            deletions.select("CHR", "del_pos", "cover_end"),
            left_on="_before",
            right_on="del_pos",
            by="CHR",
            strategy="backward",
        )
        .with_columns((pl.col("cover_end") >= pl.col("Pos")).fill_null(False).alias("under_deletion"))
        .drop("_before", "del_pos", "cover_end")
    )


def variant_info_check(lone: pl.DataFrame, vi_path) -> pl.DataFrame:  # noqa: ANN001
    vi = pl.scan_parquet(vi_path).select("CHR", "POS", "A1", "ImpFreqA1")
    keys = lone.select(pl.col("CHR").cast(pl.Int8), pl.col("Pos").alias("POS")).unique()
    at_pos = (
        vi.join(keys.lazy(), on=["CHR", "POS"], how="semi")
        .group_by("CHR", "POS")
        .agg(
            pl.col("A1").alias("vi_alleles"),
            pl.col("ImpFreqA1").alias("vi_freqs"),
        )
        .collect()
    )
    return lone.with_columns(pl.col("CHR").cast(pl.Int8)).join(
        at_pos, left_on=["CHR", "Pos"], right_on=["CHR", "POS"], how="left"
    ).with_columns(
        pl.col("vi_alleles").list.len().fill_null(0).alias("n_vi_rows"),
        pl.col("vi_alleles").list.n_unique().fill_null(0).alias("n_vi_alleles"),
        (
            pl.col("vi_alleles").list.set_difference(pl.concat_list("EA", "OA")).list.len()
            == 0
        ).alias("vi_only_own_alleles"),
        (pl.col("vi_freqs").list.sum() - 1).abs().alias("freq_sum_gap"),
    ).with_columns(
        pl.when(pl.col("n_vi_rows") == 0)
        .then(pl.lit("absent from variant_info"))
        .when(
            (pl.col("n_vi_rows") == 2)
            & (pl.col("n_vi_alleles") == 2)
            & pl.col("vi_only_own_alleles")
            & (pl.col("freq_sum_gap") <= FREQ_TOLERANCE)
        )
        .then(pl.lit("biallelic"))
        .when(pl.col("vi_only_own_alleles") & (pl.col("n_vi_rows") > pl.col("n_vi_alleles")))
        .then(pl.lit("own alleles, repeated"))
        .when(~pl.col("vi_only_own_alleles"))
        .then(pl.lit("extra alleles"))
        .otherwise(pl.lit("own alleles, freq gap"))
        .alias("vi_status")
    )


def main() -> None:
    assets = DEFAULT_RUNNER.run(
        [DECODE_SMP_AP4M1_RAW_FILE, UCSC_HG38_INDEXED_FASTA, DECODE_VARIANT_INFO]
    )
    raw = assets[DECODE_SMP_AP4M1_RAW_FILE.asset_id]
    fasta_dir = assets[UCSC_HG38_INDEXED_FASTA.asset_id]
    vi = assets[DECODE_VARIANT_INFO.asset_id]
    assert isinstance(raw, FileAsset) and isinstance(vi, FileAsset)
    assert isinstance(fasta_dir, DirectoryAsset)
    rows = add_reference_matches(read_rows(raw.path), IndexedFasta.open(fasta_dir.path))
    allele_coding(rows)

    rows = deletion_cover(rows)
    lone = rows.filter((pl.col("kind") == "snv") & (pl.col("n_at_pos") == 1)).with_columns(
        (pl.col("ImpMAF") >= 0.01).fill_null(False).alias("common")
    )
    print(f"\nlone SNVs: {lone.height} ({int(lone['common'].sum())} with MAF >= 1%)")
    print("\n(b) lone SNVs inside a deletion's span")
    print(
        lone.group_by("common")
        .agg(pl.len(), pl.col("under_deletion").sum(), pl.col("under_deletion").mean().alias("share"))
        .sort("common")
    )

    checked = variant_info_check(lone, vi.path)
    print("\n(c) variant_info at lone-SNV positions")
    print(
        checked.group_by("common", "vi_status")
        .agg(pl.len(), pl.col("under_deletion").sum().alias("under_deletion"))
        .sort("common", "len", descending=[False, True])
    )
    for status in ("extra alleles", "own alleles, freq gap", "own alleles, repeated"):
        sample = checked.filter(pl.col("vi_status") == status, pl.col("common"))
        print(f"\nexamples, {status} (MAF >= 1%): {sample.height}")
        print(
            sample.head(6).select(
                "CHR", "Pos", "EA", "OA", "ImpMAF", "under_deletion", "vi_alleles", "vi_freqs"
            )
        )


if __name__ == "__main__":
    main()
