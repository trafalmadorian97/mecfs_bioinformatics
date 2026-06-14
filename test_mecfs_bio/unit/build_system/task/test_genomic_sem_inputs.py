"""
Unit tests for the rpy2-free input helpers in ``_genomic_sem_inputs`` used by
the full-Python GenomicSEM path: sample size / prevalence extraction and the
gwaslab -> canonical-munge dataframe conversion (``build_munge_input_df``).
"""

from pathlib import Path

import polars as pl
import polars.testing
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
    QuantPhenotype,
)
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
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_inputs import (
    add_sample_size_if_missing,
    build_munge_input_df,
    get_prevs,
    get_sample_size,
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


def _make_dummy_sumstats(n_rows: int = 5) -> pl.DataFrame:
    return pl.DataFrame(
        {
            GWASLAB_RSID_COL: [f"rs{i}" for i in range(n_rows)],
            GWASLAB_EFFECT_ALLELE_COL: ["A"] * n_rows,
            GWASLAB_NON_EFFECT_ALLELE_COL: ["G"] * n_rows,
            GWASLAB_BETA_COL: [0.1 * i for i in range(n_rows)],
            GWASLAB_SE_COL: [0.05] * n_rows,
            GWASLAB_P_COL: [0.5] * n_rows,
            GWASLAB_SAMPLE_SIZE_COLUMN: [10000] * n_rows,
            GWASLAB_EFFECT_ALLELE_FREQ_COL: [0.3] * n_rows,
        }
    )


def _make_source(
    asset_id: str,
    alias: str,
    sample_info,
) -> GenomicSEMSumstatsSource:
    task = FakeTask(
        SimpleFileMeta(
            AssetId(asset_id),
            read_spec=DataFrameReadSpec(DataFrameTextFormat(",")),
        ),
    )
    return GenomicSEMSumstatsSource(task=task, alias=alias, sample_info=sample_info)


def _csv_source(df: pl.DataFrame, tmp_path: Path, alias: str = "trait"):
    """Write df as CSV and return (source, fetch) reading it back."""
    source_path = tmp_path / f"{alias}.csv"
    df.write_csv(source_path)
    source = _make_source(alias, alias, QuantPhenotype(total_sample_size=10000))

    def fetch(asset_id: AssetId) -> Asset:
        return FileAsset(source_path)

    return source, fetch


def test_get_sample_size_quant():
    info = QuantPhenotype(total_sample_size=12345)
    src = _make_source("x", "x", info)
    assert get_sample_size(src) == pytest.approx(12345.0)


def test_get_sample_size_binary():
    info = BinaryPhenotypeSampleInfo(
        sample_prevalence=0.3,
        estimated_population_prevalence=0.05,
        total_sample_size=4000,
    )
    src = _make_source("x", "x", info)
    assert get_sample_size(src) == pytest.approx(4000.0)


def test_get_prevs_quant_returns_nan():
    samp, pop = get_prevs(QuantPhenotype(total_sample_size=10))
    assert samp != samp  # NaN
    assert pop != pop


def test_get_prevs_binary():
    samp, pop = get_prevs(
        BinaryPhenotypeSampleInfo(
            sample_prevalence=0.2,
            estimated_population_prevalence=0.01,
            total_sample_size=5000,
        )
    )
    assert samp == pytest.approx(0.2)
    assert pop == pytest.approx(0.01)


def test_add_sample_size_if_missing_uses_quant_total():
    df = _make_dummy_sumstats(n_rows=3).drop(GWASLAB_SAMPLE_SIZE_COLUMN)
    out = add_sample_size_if_missing(
        df, sample_info=QuantPhenotype(total_sample_size=999)
    )
    assert (out[GWASLAB_SAMPLE_SIZE_COLUMN] == 999).all()


def test_add_sample_size_if_missing_keeps_existing_column():
    df = _make_dummy_sumstats(n_rows=3)
    out = add_sample_size_if_missing(
        df, sample_info=QuantPhenotype(total_sample_size=1)
    )
    pl.testing.assert_frame_equal(out, df)


def test_add_sample_size_if_missing_raises_when_missing_total():
    df = _make_dummy_sumstats(n_rows=3).drop(GWASLAB_SAMPLE_SIZE_COLUMN)
    with pytest.raises(ValueError):
        add_sample_size_if_missing(
            df, sample_info=QuantPhenotype(total_sample_size=None)
        )


def test_build_munge_input_df_renames_to_canonical_columns(tmp_path: Path):
    """gwaslab columns are mapped to the canonical munge names munge/sumstats use."""
    source, fetch = _csv_source(_make_dummy_sumstats(n_rows=4), tmp_path)

    out = build_munge_input_df(source, fetch)

    assert set(out.columns) == {
        MUNGE_SNP_COL,
        MUNGE_A1_COL,
        MUNGE_A2_COL,
        MUNGE_EFFECT_COL,
        MUNGE_SE_COL,
        MUNGE_P_COL,
        MUNGE_N_COL,
        MUNGE_MAF_COL,
    }
    assert len(out) == 4
    assert out[MUNGE_SNP_COL].to_list() == [f"rs{i}" for i in range(4)]
    assert (out[MUNGE_N_COL] == 10000).all()


def test_build_munge_input_df_omits_maf_when_freq_missing(tmp_path: Path):
    """EAF is optional; absent it, MAF is omitted rather than filled with NaN."""
    df = _make_dummy_sumstats(n_rows=3).drop(GWASLAB_EFFECT_ALLELE_FREQ_COL)
    source, fetch = _csv_source(df, tmp_path)

    out = build_munge_input_df(source, fetch)
    assert MUNGE_MAF_COL not in out.columns


def test_build_munge_input_df_raises_when_required_column_missing(tmp_path: Path):
    df = _make_dummy_sumstats(n_rows=3).drop(GWASLAB_BETA_COL)
    source, fetch = _csv_source(df, tmp_path)

    with pytest.raises(AssertionError):
        build_munge_input_df(source, fetch)


def test_build_munge_input_df_casts_categorical_alleles_to_string(tmp_path: Path):
    """
    Parquet sources deliver SNP/allele columns as Categorical (dictionary-
    encoded), which breaks the `.str` ops and SNP joins in munge/sumstats.
    build_munge_input_df must normalise SNP/A1/A2 to String regardless.
    """
    df = _make_dummy_sumstats(n_rows=4).with_columns(
        pl.col(GWASLAB_RSID_COL).cast(pl.Categorical),
        pl.col(GWASLAB_EFFECT_ALLELE_COL).cast(pl.Categorical),
        pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).cast(pl.Categorical),
    )
    source_path = tmp_path / "source.parquet"
    df.write_parquet(source_path)

    task = FakeTask(
        SimpleFileMeta(
            AssetId("trait_cat"),
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        ),
    )
    source = GenomicSEMSumstatsSource(
        task=task,
        alias="trait_cat",
        sample_info=QuantPhenotype(total_sample_size=10000),
    )

    def fetch(asset_id: AssetId) -> Asset:
        return FileAsset(source_path)

    out = build_munge_input_df(source, fetch)
    assert out.schema[MUNGE_SNP_COL] == pl.String
    assert out.schema[MUNGE_A1_COL] == pl.String
    assert out.schema[MUNGE_A2_COL] == pl.String
