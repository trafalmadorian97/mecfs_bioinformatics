"""
rpy2-free input-preparation helpers used by the full-Python GenomicSEM tasks:
deriving the sample size, reading a fetched tabular asset, building the canonical
munge-format DataFrame, and resolving the LD reference directory.

R-only input helpers (the munge-TSV writer, single-file path resolution, source
validation, sumstats method flags, lavaan-component sanitisation) live in
``_genomic_sem_r_inputs``.
"""

from __future__ import annotations

import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import polars as pl

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_config import (
    MUNGE_A1_COL,
    MUNGE_A2_COL,
    MUNGE_EFFECT_COL,
    MUNGE_MAF_COL,
    MUNGE_N_COL,
    MUNGE_P_COL,
    MUNGE_SE_COL,
    MUNGE_SNP_COL,
    GenomicSEMSumstatsSource,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_P_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)


def get_sample_size(source: GenomicSEMSumstatsSource) -> float:
    """Total sample size for the source, or NaN ("not provided", so the file's
    own N column stands) when the source carries no ``sample_size``."""
    if source.sample_size is not None:
        return float(source.sample_size)
    return float("nan")


def read_dataframe_from_task(task: Task, fetch: Fetch) -> pl.DataFrame:
    """
    Read a fetched single-file tabular asset as a polars DataFrame using the
    format declared in the task's metadata read_spec. This lets the upstream
    task choose whatever on-disk format is convenient (tsv, space-delimited,
    parquet, ...) without the reader needing to know which.
    """
    asset = fetch(task.asset_id)
    return scan_dataframe_asset(asset, task.meta).collect().to_polars()


def build_munge_input_df(
    source: GenomicSEMSumstatsSource, fetch: Fetch
) -> pl.DataFrame:
    """
    Read the source sumstats (via its read_spec) and rename to the canonical
    munge columns (SNP, A1, A2, effect, SE, P, N, and MAF when present). This
    is what ``munge_sumstats`` / ``sumstats`` consume; ``write_munge_input``
    additionally serialises it to disk for the R workflow.
    """
    asset = fetch(source.asset_id)
    df = (
        source.pipe.process(scan_dataframe_asset(asset, source.task.meta))
        .collect()
        .to_polars()
    )
    df = add_sample_size_if_missing(df, sample_size=source.sample_size)

    required = [
        GWASLAB_RSID_COL,
        GWASLAB_EFFECT_ALLELE_COL,
        GWASLAB_NON_EFFECT_ALLELE_COL,
        GWASLAB_BETA_COL,
        GWASLAB_SE_COL,
        GWASLAB_P_COL,
        GWASLAB_SAMPLE_SIZE_COLUMN,
    ]
    missing = [c for c in required if c not in df.columns]
    assert not missing, f"Source '{source.alias}' is missing required columns {missing}"

    # Cast the SNP / allele columns to String. Parquet sources can sometimes deliver them as
    # Categorical (dictionary-encoded), which breaks the `.str` operations and
    # SNP joins in munge/sumstats; normalising here keeps those ports simple.
    select_exprs: list[pl.Expr] = [
        pl.col(GWASLAB_RSID_COL).cast(pl.String).alias(MUNGE_SNP_COL),
        pl.col(GWASLAB_EFFECT_ALLELE_COL).cast(pl.String).alias(MUNGE_A1_COL),
        pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).cast(pl.String).alias(MUNGE_A2_COL),
        pl.col(GWASLAB_BETA_COL).alias(MUNGE_EFFECT_COL),
        pl.col(GWASLAB_SE_COL).alias(MUNGE_SE_COL),
        pl.col(GWASLAB_P_COL).alias(MUNGE_P_COL),
        pl.col(GWASLAB_SAMPLE_SIZE_COLUMN).alias(MUNGE_N_COL),
    ]
    if GWASLAB_EFFECT_ALLELE_FREQ_COL in df.columns:
        # Fold the effect-allele frequency to the minor-allele frequency so the
        # MUNGE_MAF_COL column actually holds a MAF. The source reports the
        # frequency of the effect allele, which need not be the minor allele;
        # downstream munge/sumstats only ever use the minor frequency (for the
        # MAF filter and varSNP = 2p(1-p)), so we normalise it at the source.
        select_exprs.append(
            pl.when(pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL) <= 0.5)
            .then(pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL))
            .otherwise(1.0 - pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL))
            .alias(MUNGE_MAF_COL)
        )
    return df.select(select_exprs)


def add_sample_size_if_missing(
    df: pl.DataFrame, sample_size: float | None
) -> pl.DataFrame:
    if GWASLAB_SAMPLE_SIZE_COLUMN in df.columns:
        return df
    if sample_size is None:
        raise ValueError(
            "sample_size must be set when the source dataframe does not include "
            "an N column"
        )
    return df.with_columns(pl.lit(sample_size).alias(GWASLAB_SAMPLE_SIZE_COLUMN))


@contextmanager
def ld_dir_with_genomic_sem_naming(
    ld_path: Path, basename_prefix: str
) -> Iterator[Path]:
    """
    GenomicSEM::ldsc constructs LD score paths as "<ld>/<chr>.l2.ldscore.gz",
    so when the on-disk files have a basename prefix (e.g. "LDscore.") we
    build a tempdir of symlinks with the prefix stripped and yield that.
    """
    if not basename_prefix:
        yield ld_path
        return
    with tempfile.TemporaryDirectory() as tmp_str:
        tmp = Path(tmp_str)
        for src in ld_path.iterdir():
            if not src.name.startswith(basename_prefix):
                continue
            stripped = src.name[len(basename_prefix) :]
            (tmp / stripped).symlink_to(src.resolve())
        yield tmp


def resolve_ld_path(ld_ref_task: Task, fetch: Fetch) -> Path:
    asset = fetch(ld_ref_task.asset_id)
    assert isinstance(asset, DirectoryAsset)
    # Absolute path so the value remains correct after R chdirs during munge.
    # The basename-prefix in munge_config.ld_file_filename_pattern is applied
    # downstream via a symlink farm — do not concatenate it here.
    return asset.path.resolve()
