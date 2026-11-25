"""
How far do deCODE (Iceland, variant_info) and DecodeME (UK, harmonized GRCh38) effect allele
frequencies differ at the same variant? Sets the threshold of a frequency-concordance check in
coloc alignment.

deCODE side: every clean row of the AP4M1 reference file (keep_clean_rows), EAF resolved from
the clean variant_info. deCODE's otherAllele is the GRCh38 reference, as is the harmonized
DecodeME NEA, so the two are joined exactly on (CHR, POS, EA, NEA). Reports the distribution of
d = EAF_decode - EAF_decodeme per kind (SNV / indel, palindromic or not) and MAF band, the tail
counts, and how many tail rows would instead agree with the flipped frequency (1 - EAF).

  pixi r python experiments/claude/decode_proteomics_survey/eaf_concordance_decode_vs_decodeme.py
"""

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.decode_somascan.raw.decode_reference_protein_file import (
    DECODE_SMP_AP4M1_RAW_FILE,
)
from mecfs_bio.assets.gwas.decode_somascan.raw.decode_variant_info import (
    DECODE_CLEAN_VARIANT_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_build_38_harmonized import (
    DECODE_ME_GWAS_1_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.task.decode_database.build_decode_variant_info_task import (
    keep_clean_rows,
)
from mecfs_bio.build_system.task.decode_database.decode_sumstats_format import (
    chrom_code_expr,
)
from mecfs_bio.build_system.task.decode_database.effect_allele_frequency import (
    resolve_effect_allele_frequency,
)
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_rsid import (
    is_palindromic_expr,
)

pl.Config.set_tbl_rows(40)
pl.Config.set_tbl_cols(-1)
pl.Config.set_tbl_width_chars(220)
KEY = ["CHR", "POS", "EA", "NEA"]
TAILS = (0.05, 0.1, 0.15, 0.2, 0.3, 0.4)


def decode_rows(raw_path, clean_path) -> pl.DataFrame:  # noqa: ANN001
    rows = pl.read_csv(
        raw_path,
        separator="\t",
        columns=["Chrom", "Pos", "effectAllele", "otherAllele", "ImpMAF"],
        null_values=["."],
        quote_char=None,
        schema_overrides={"Pos": pl.Int64, "ImpMAF": pl.Float64},
    ).select(
        chrom_code_expr("Chrom").cast(pl.Int8).alias("CHR"),
        pl.col("Pos").alias("POS"),
        pl.col("effectAllele").alias("EA"),
        pl.col("otherAllele").alias("NEA"),
        "ImpMAF",
    )
    clean_vi = pl.read_parquet(clean_path)
    clean = keep_clean_rows(rows, clean_vi)
    print(f"deCODE rows {rows.height}, clean {clean.height}", flush=True)
    eaf = resolve_effect_allele_frequency(clean, clean_vi)
    return clean.with_columns(eaf.alias("EAF_decode")).drop_nulls("EAF_decode")


def main() -> None:
    assets = DEFAULT_RUNNER.run(
        [DECODE_SMP_AP4M1_RAW_FILE, DECODE_CLEAN_VARIANT_INFO, DECODE_ME_GWAS_1_BUILD_38_HARMONIZED]
    )
    raw = assets[DECODE_SMP_AP4M1_RAW_FILE.asset_id]
    clean = assets[DECODE_CLEAN_VARIANT_INFO.asset_id]
    assert isinstance(raw, FileAsset) and isinstance(clean, FileAsset)
    decode = decode_rows(raw.path, clean.path)
    trait = (
        scan_dataframe_asset(
            assets[DECODE_ME_GWAS_1_BUILD_38_HARMONIZED.asset_id],
            meta=DECODE_ME_GWAS_1_BUILD_38_HARMONIZED.meta,
            parquet_backend="polars",
        )
        .to_native()
        .select(
            pl.col("CHR").cast(pl.Int8),
            pl.col("POS").cast(pl.Int64),
            pl.col("EA").cast(pl.String),
            pl.col("NEA").cast(pl.String),
            pl.col("EAF").cast(pl.Float64).alias("EAF_trait"),
        )
        .collect()
    )
    print(f"DecodeME rows {trait.height}, EAF null {trait['EAF_trait'].null_count()}")
    d = pl.col("EAF_decode") - pl.col("EAF_trait")
    joined = (
        decode.join(trait.drop_nulls("EAF_trait"), on=KEY, how="inner")
        .with_columns(
            d.alias("d"),
            d.abs().alias("abs_d"),
            (pl.col("EAF_decode") - (1 - pl.col("EAF_trait"))).abs().alias("abs_d_flipped"),
            pl.when((pl.col("EA").str.len_bytes() == 1) & (pl.col("NEA").str.len_bytes() == 1))
            .then(
                pl.when(is_palindromic_expr("EA", "NEA"))
                .then(pl.lit("snv_palindromic"))
                .otherwise(pl.lit("snv"))
            )
            .otherwise(pl.lit("indel"))
            .alias("kind"),
            pl.min_horizontal("EAF_trait", 1 - pl.col("EAF_trait"))
            .cut([0.01, 0.05, 0.2])
            .alias("maf_trait"),
        )
    )
    print(f"matched {joined.height}; corr {joined.select(pl.corr('EAF_decode', 'EAF_trait')).item():.4f}")
    q = [0.5, 0.9, 0.99, 0.999, 0.9999]
    print(
        joined.group_by("kind", "maf_trait")
        .agg(
            pl.len().alias("n"),
            *[pl.col("abs_d").quantile(p).alias(f"q{p}") for p in q],
            pl.col("abs_d").max().alias("max"),
        )
        .sort("kind", "maf_trait")
    )
    print(
        joined.group_by("kind")
        .agg(
            pl.len().alias("n"),
            *[(pl.col("abs_d") > t).sum().alias(f">{t}") for t in TAILS],
            *[
                ((pl.col("abs_d") > t) & (pl.col("abs_d_flipped") < 0.05)).sum().alias(f">{t}&flip<.05")
                for t in (0.2, 0.3)
            ],
        )
        .sort("kind")
    )
    print(
        joined.filter(pl.col("abs_d") > 0.3)
        .sort("abs_d", descending=True)
        .select(*KEY, "ImpMAF", "EAF_decode", "EAF_trait", "abs_d", "abs_d_flipped", "kind")
        .head(25)
    )


if __name__ == "__main__":
    main()
