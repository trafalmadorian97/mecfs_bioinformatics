"""
rpy2-free input-preparation helpers shared across the GenomicSEM task family:
reading/writing the munge-format inputs, deriving sample sizes and
prevalences, resolving fetched asset paths, and validating sources.

These helpers contain no R calls, so both the R-backed tasks and the
full-Python tasks can rely on them.
"""

from __future__ import annotations

import re
import tempfile
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path

import polars as pl
import structlog

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
    PhenotypeInfo,
    QuantPhenotype,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_config import (
    LINEAR_PROB,
    LOGISTIC,
    MUNGE_A1_COL,
    MUNGE_A2_COL,
    MUNGE_EFFECT_COL,
    MUNGE_MAF_COL,
    MUNGE_N_COL,
    MUNGE_P_COL,
    MUNGE_SE_COL,
    MUNGE_SNP_COL,
    OLS,
    GenomicSEMConfig,
    GenomicSEMGWASSumstatsSource,
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

logger = structlog.get_logger()


def _get_sample_size(source: GenomicSEMSumstatsSource) -> float:
    info = source.sample_info
    if isinstance(info, BinaryPhenotypeSampleInfo):
        if info.total_sample_size is not None:
            return float(info.total_sample_size)
        return float("nan")
    if isinstance(info, QuantPhenotype):
        if info.total_sample_size is not None:
            return float(info.total_sample_size)
        return float("nan")
    raise ValueError(f"Unknown sample info type {type(info)}")


def _get_prevs(info: PhenotypeInfo) -> tuple[float, float]:
    if isinstance(info, BinaryPhenotypeSampleInfo):
        return info.sample_prevalence, info.estimated_population_prevalence
    if isinstance(info, QuantPhenotype):
        return float("nan"), float("nan")
    raise ValueError(f"Unknown sample info type {type(info)}")


def _write_munge_input(
    source: GenomicSEMSumstatsSource, fetch: Fetch, tmp_dir: Path
) -> Path:
    """
    Read the source sumstats and write a TSV with the column names munge expects.
    """
    asset = fetch(source.asset_id)
    df = (
        source.pipe.process(scan_dataframe_asset(asset, source.task.meta))
        .collect()
        .to_polars()
    )
    df = _add_sample_size_if_missing(df, sample_info=source.sample_info)

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

    select_exprs: list[pl.Expr] = [
        pl.col(GWASLAB_RSID_COL).alias(MUNGE_SNP_COL),
        pl.col(GWASLAB_EFFECT_ALLELE_COL).alias(MUNGE_A1_COL),
        pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).alias(MUNGE_A2_COL),
        pl.col(GWASLAB_BETA_COL).alias(MUNGE_EFFECT_COL),
        pl.col(GWASLAB_SE_COL).alias(MUNGE_SE_COL),
        pl.col(GWASLAB_P_COL).alias(MUNGE_P_COL),
        pl.col(GWASLAB_SAMPLE_SIZE_COLUMN).alias(MUNGE_N_COL),
    ]
    if GWASLAB_EFFECT_ALLELE_FREQ_COL in df.columns:
        select_exprs.append(pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL).alias(MUNGE_MAF_COL))
    munge_df = df.select(select_exprs)

    output_path = tmp_dir / f"{source.alias}.sumstats.txt"
    munge_df.write_csv(output_path, separator="\t")
    logger.debug(
        f"Wrote munge input for '{source.alias}' to {output_path} "
        f"({len(munge_df)} variants)"
    )
    return output_path


def _add_sample_size_if_missing(
    df: pl.DataFrame, sample_info: PhenotypeInfo
) -> pl.DataFrame:
    if GWASLAB_SAMPLE_SIZE_COLUMN in df.columns:
        return df
    if isinstance(sample_info, BinaryPhenotypeSampleInfo):
        if sample_info.total_sample_size is None:
            raise ValueError(
                "BinaryPhenotypeSampleInfo.total_sample_size must be set when "
                "the source dataframe does not include an N column"
            )
        n = sample_info.total_sample_size
    elif isinstance(sample_info, QuantPhenotype):
        if sample_info.total_sample_size is None:
            raise ValueError(
                "QuantPhenotype.total_sample_size must be set when the source "
                "dataframe does not include an N column"
            )
        n = sample_info.total_sample_size
    else:
        raise ValueError(f"Unknown sample info type {type(sample_info)}")
    return df.with_columns(pl.lit(n).alias(GWASLAB_SAMPLE_SIZE_COLUMN))


@contextmanager
def _ld_dir_with_genomic_sem_naming(
    ld_path: str, basename_prefix: str
) -> Iterator[str]:
    """
    GenomicSEM::ldsc constructs LD score paths as "<ld>/<chr>.l2.ldscore.gz",
    so when the on-disk files have a basename prefix (e.g. "LDscore.") we
    build a tempdir of symlinks with the prefix stripped and yield that.
    """
    if not basename_prefix:
        yield ld_path
        return
    src_dir = Path(ld_path)
    with tempfile.TemporaryDirectory() as tmp_str:
        tmp = Path(tmp_str)
        for src in src_dir.iterdir():
            if not src.name.startswith(basename_prefix):
                continue
            stripped = src.name[len(basename_prefix) :]
            (tmp / stripped).symlink_to(src.resolve())
        yield str(tmp)


def _resolve_ld_path(
    ld_ref_task: Task, fetch: Fetch, munge_config: GenomicSEMConfig
) -> str:
    asset = fetch(ld_ref_task.asset_id)
    assert isinstance(asset, DirectoryAsset)
    # Absolute path so the value remains correct after R chdirs during munge.
    # The basename-prefix in munge_config.ld_file_filename_pattern is applied
    # downstream via a symlink farm — do not concatenate it here.
    return str(asset.path.resolve())


def _resolve_file_path(task: Task, fetch: Fetch) -> str:
    asset = fetch(task.asset_id)
    assert isinstance(asset, FileAsset)
    return str(asset.path.resolve())


def _validate_sources(sources: Sequence[GenomicSEMGWASSumstatsSource]) -> None:
    assert len(sources) >= 2, "GenomicSEM GWAS requires at least two traits"
    aliases = [s.alias for s in sources]
    assert len(set(aliases)) == len(aliases), (
        f"Trait aliases must be unique. Got: {aliases}"
    )


def _gwas_method_flags(
    sources: Sequence[GenomicSEMGWASSumstatsSource],
) -> tuple[list[bool], list[bool], list[bool]]:
    """
    GenomicSEM::sumstats expects a (se.logit, OLS, linprob) BoolVector triple.
    Exactly one is TRUE per trait.
    """
    se_logit: list[bool] = []
    ols: list[bool] = []
    linprob: list[bool] = []
    for s in sources:
        se_logit.append(s.gwas_method == LOGISTIC)
        ols.append(s.gwas_method == OLS)
        linprob.append(s.gwas_method == LINEAR_PROB)
    return se_logit, ols, linprob


def _sanitize_component_name(name: str) -> str:
    """Map a lavaan sub-component (e.g. 'F1~SNP') to a filename stem ('F1_SNP')."""
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")
