"""
Unit tests for GenomicSEMTask.

These tests cover the Python-side wiring of the task: dependency graph,
metadata, source-file munge formatting, sample-size handling, and prevalence
extraction. Running the actual GenomicSEM R workflow against a real LD
reference is left to a future system test.
"""

from pathlib import Path

import polars as pl
import polars.testing
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
    QuantPhenotype,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_task import (
    MUNGE_A1_COL,
    MUNGE_A2_COL,
    MUNGE_EFFECT_COL,
    MUNGE_MAF_COL,
    MUNGE_N_COL,
    MUNGE_P_COL,
    MUNGE_SE_COL,
    MUNGE_SNP_COL,
    GenomicSEMConfig,
    GenomicSEMSumstatsSource,
    GenomicSEMTask,
    _add_sample_size_if_missing,
    _get_prevs,
    _get_sample_size,
    _write_munge_input,
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


def test_create_meta_and_deps():
    """
    The factory should construct meta describing a multi-trait result directory
    and deps should include the LD ref task, hapmap task, and each source task.
    """
    ld_ref_task = FakeTask.create_with_filemeta("ld_ref")
    hapmap_task = FakeTask.create_with_filemeta("hapmap")
    src_a = _make_source("trait_a_data", "a", QuantPhenotype(total_sample_size=10000))
    src_b = _make_source("trait_b_data", "b", QuantPhenotype(total_sample_size=20000))

    task = GenomicSEMTask.create(
        asset_id="genomic_sem_test",
        sources=[src_a, src_b],
        ld_ref_task=ld_ref_task,
        hapmap_snps_task=hapmap_task,
        factor_model="F1 =~ a + b",
    )

    assert isinstance(task.meta, ResultDirectoryMeta)
    assert task.meta.asset_id == "genomic_sem_test"
    assert task.meta.project == "genomic_sem"
    assert task.deps == [ld_ref_task, hapmap_task, src_a.task, src_b.task]


def test_requires_at_least_two_sources():
    """
    GenomicSEM does multi-trait analysis; passing fewer than two sources is a bug.
    """
    src_a = _make_source("trait_a_data", "a", QuantPhenotype(total_sample_size=10000))
    with pytest.raises(AssertionError):
        GenomicSEMTask.create(
            asset_id="genomic_sem_test",
            sources=[src_a],
            ld_ref_task=FakeTask.create_with_filemeta("ld_ref"),
            hapmap_snps_task=FakeTask.create_with_filemeta("hapmap"),
            factor_model="F1 =~ a",
        )


def test_rejects_duplicate_aliases():
    """
    Trait aliases double as munge output filenames so they must be unique.
    """
    src_a = _make_source("trait_a_data", "a", QuantPhenotype(total_sample_size=10000))
    src_b = _make_source("trait_b_data", "a", QuantPhenotype(total_sample_size=20000))
    with pytest.raises(AssertionError):
        GenomicSEMTask.create(
            asset_id="genomic_sem_test",
            sources=[src_a, src_b],
            ld_ref_task=FakeTask.create_with_filemeta("ld_ref"),
            hapmap_snps_task=FakeTask.create_with_filemeta("hapmap"),
            factor_model="F1 =~ a",
        )


def test_get_sample_size_quant():
    info = QuantPhenotype(total_sample_size=12345)
    src = _make_source("x", "x", info)
    assert _get_sample_size(src) == pytest.approx(12345.0)


def test_get_sample_size_binary():
    info = BinaryPhenotypeSampleInfo(
        sample_prevalence=0.3,
        estimated_population_prevalence=0.05,
        total_sample_size=4000,
    )
    src = _make_source("x", "x", info)
    assert _get_sample_size(src) == pytest.approx(4000.0)


def test_get_prevs_quant_returns_nan():
    samp, pop = _get_prevs(QuantPhenotype(total_sample_size=10))
    assert samp != samp  # NaN
    assert pop != pop


def test_get_prevs_binary():
    samp, pop = _get_prevs(
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
    out = _add_sample_size_if_missing(
        df, sample_info=QuantPhenotype(total_sample_size=999)
    )
    assert (out[GWASLAB_SAMPLE_SIZE_COLUMN] == 999).all()


def test_add_sample_size_if_missing_keeps_existing_column():
    df = _make_dummy_sumstats(n_rows=3)
    out = _add_sample_size_if_missing(
        df, sample_info=QuantPhenotype(total_sample_size=1)
    )
    pl.testing.assert_frame_equal(out, df)


def test_add_sample_size_if_missing_raises_when_missing_total():
    df = _make_dummy_sumstats(n_rows=3).drop(GWASLAB_SAMPLE_SIZE_COLUMN)
    with pytest.raises(ValueError):
        _add_sample_size_if_missing(
            df, sample_info=QuantPhenotype(total_sample_size=None)
        )


def test_write_munge_input_renames_columns_and_writes_tsv(tmp_path: Path):
    """
    _write_munge_input is the bridge between gwaslab format and what
    GenomicSEM::munge expects. Verify the file lands on disk with the right
    columns.
    """
    df = _make_dummy_sumstats(n_rows=4)
    source_path = tmp_path / "source.csv"
    df.write_csv(source_path)

    task = FakeTask(
        SimpleFileMeta(
            AssetId("trait_a"),
            read_spec=DataFrameReadSpec(DataFrameTextFormat(",")),
        ),
    )
    source = GenomicSEMSumstatsSource(
        task=task,
        alias="trait_a",
        sample_info=QuantPhenotype(total_sample_size=10000),
    )

    def fetch(asset_id: AssetId) -> Asset:
        assert asset_id == "trait_a"
        return FileAsset(source_path)

    output_path = _write_munge_input(source=source, fetch=fetch, tmp_dir=tmp_path)

    assert output_path == tmp_path / "trait_a.sumstats.txt"
    assert output_path.is_file()
    written = pl.read_csv(output_path, separator="\t")
    assert set(written.columns) == {
        MUNGE_SNP_COL,
        MUNGE_A1_COL,
        MUNGE_A2_COL,
        MUNGE_EFFECT_COL,
        MUNGE_SE_COL,
        MUNGE_P_COL,
        MUNGE_N_COL,
        MUNGE_MAF_COL,
    }
    assert len(written) == 4
    assert written[MUNGE_SNP_COL].to_list() == [f"rs{i}" for i in range(4)]
    assert (written[MUNGE_N_COL] == 10000).all()


def test_write_munge_input_omits_maf_when_freq_missing(tmp_path: Path):
    """
    EAF is optional; if it's not in the source dataframe, MAF should be absent
    from the munge input rather than filled with NaN.
    """
    df = _make_dummy_sumstats(n_rows=3).drop(GWASLAB_EFFECT_ALLELE_FREQ_COL)
    source_path = tmp_path / "source.csv"
    df.write_csv(source_path)

    task = FakeTask(
        SimpleFileMeta(
            AssetId("trait_x"),
            read_spec=DataFrameReadSpec(DataFrameTextFormat(",")),
        ),
    )
    source = GenomicSEMSumstatsSource(
        task=task,
        alias="trait_x",
        sample_info=QuantPhenotype(total_sample_size=10000),
    )

    def fetch(asset_id: AssetId) -> Asset:
        return FileAsset(source_path)

    output_path = _write_munge_input(source=source, fetch=fetch, tmp_dir=tmp_path)
    written = pl.read_csv(output_path, separator="\t")
    assert MUNGE_MAF_COL not in written.columns


def test_write_munge_input_raises_when_required_column_missing(tmp_path: Path):
    df = _make_dummy_sumstats(n_rows=3).drop(GWASLAB_BETA_COL)
    source_path = tmp_path / "source.csv"
    df.write_csv(source_path)

    task = FakeTask(
        SimpleFileMeta(
            AssetId("trait_x"),
            read_spec=DataFrameReadSpec(DataFrameTextFormat(",")),
        ),
    )
    source = GenomicSEMSumstatsSource(
        task=task,
        alias="trait_x",
        sample_info=QuantPhenotype(total_sample_size=10000),
    )

    def fetch(asset_id: AssetId) -> Asset:
        return FileAsset(source_path)

    with pytest.raises(AssertionError):
        _write_munge_input(source=source, fetch=fetch, tmp_dir=tmp_path)


def test_config_has_sensible_defaults():
    config = GenomicSEMConfig()
    assert config.estimation == "DWLS"
    assert config.info_filter == pytest.approx(0.9)
    assert config.maf_filter == pytest.approx(0.01)


def test_task_uses_simple_directory_meta_when_passed(tmp_path: Path):
    """
    Constructing GenomicSEMTask directly (rather than via .create) should
    accept any Meta — exercise the SimpleDirectoryMeta path used by tests.
    """
    src_a = _make_source("a_data", "a", QuantPhenotype(total_sample_size=10000))
    src_b = _make_source("b_data", "b", QuantPhenotype(total_sample_size=10000))
    task = GenomicSEMTask(
        meta=SimpleDirectoryMeta(AssetId("gsem_dir")),
        sources=[src_a, src_b],
        ld_ref_task=FakeTask.create_with_filemeta("ld_ref"),
        hapmap_snps_task=FakeTask.create_with_filemeta("hapmap"),
        factor_model="F1 =~ a + b",
    )
    assert task.asset_id == "gsem_dir"
