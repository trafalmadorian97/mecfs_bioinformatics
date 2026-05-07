"""
Unit tests for GenomicSEMCommonFactorGWASPythonTask -- the Python-kernel
sibling of GenomicSEMCommonFactorGWASTask. Tests cover Python-side wiring
only (deps, meta, k=2 assertion). The numerical correctness of the kernel
itself lives in test_common_factor_kernel.py.
"""

import pytest

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
    QuantPhenotype,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_common_factor_gwas_python_task import (
    GenomicSEMCommonFactorGWASPythonTask,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_task import (
    GenomicSEMSumstatsSource,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_user_gwas_task import (
    OLS,
    GenomicSEMGWASSumstatsSource,
)


def _make_gwas_source(asset_id: str, alias: str) -> GenomicSEMGWASSumstatsSource:
    inner_task = FakeTask(
        SimpleFileMeta(
            AssetId(asset_id),
            read_spec=DataFrameReadSpec(DataFrameTextFormat(",")),
        ),
    )
    inner = GenomicSEMSumstatsSource(
        task=inner_task,
        alias=alias,
        sample_info=QuantPhenotype(total_sample_size=10000),
    )
    return GenomicSEMGWASSumstatsSource(source=inner, gwas_method=OLS)


def test_create_meta_and_deps():
    ld_ref = FakeTask.create_with_filemeta("ld_ref")
    hapmap = FakeTask.create_with_filemeta("hapmap")
    ref = FakeTask.create_with_filemeta("sumstats_ref")
    src_a = _make_gwas_source("trait_a_data", "a")
    src_b = _make_gwas_source("trait_b_data", "b")

    task = GenomicSEMCommonFactorGWASPythonTask.create(
        asset_id="cf_gwas_py",
        sources=[src_a, src_b],
        ld_ref_task=ld_ref,
        hapmap_snps_task=hapmap,
        sumstats_ref_task=ref,
    )

    assert isinstance(task.meta, ResultDirectoryMeta)
    assert task.meta.asset_id == "cf_gwas_py"
    assert task.meta.project == "genomic_sem"
    assert task.deps == [ld_ref, hapmap, ref, src_a.task, src_b.task]


def test_requires_at_least_two_sources():
    src_a = _make_gwas_source("trait_a_data", "a")
    with pytest.raises(AssertionError):
        GenomicSEMCommonFactorGWASPythonTask.create(
            asset_id="cf_gwas_py",
            sources=[src_a],
            ld_ref_task=FakeTask.create_with_filemeta("ld_ref"),
            hapmap_snps_task=FakeTask.create_with_filemeta("hapmap"),
            sumstats_ref_task=FakeTask.create_with_filemeta("sumstats_ref"),
        )


def test_v1_rejects_more_than_two_traits():
    """
    v1 is k=2 only; constructing with k=3 must fail at attrs validation,
    not silently produce a bad result later.
    """
    src_a = _make_gwas_source("trait_a_data", "a")
    src_b = _make_gwas_source("trait_b_data", "b")
    src_c = _make_gwas_source("trait_c_data", "c")
    with pytest.raises(AssertionError, match="k=2"):
        GenomicSEMCommonFactorGWASPythonTask.create(
            asset_id="cf_gwas_py",
            sources=[src_a, src_b, src_c],
            ld_ref_task=FakeTask.create_with_filemeta("ld_ref"),
            hapmap_snps_task=FakeTask.create_with_filemeta("hapmap"),
            sumstats_ref_task=FakeTask.create_with_filemeta("sumstats_ref"),
        )


def test_uses_simple_directory_meta_when_passed():
    src_a = _make_gwas_source("a_data", "a")
    src_b = _make_gwas_source("b_data", "b")
    task = GenomicSEMCommonFactorGWASPythonTask(
        meta=SimpleDirectoryMeta(AssetId("cf_dir_py")),
        sources=[src_a, src_b],
        ld_ref_task=FakeTask.create_with_filemeta("ld_ref"),
        hapmap_snps_task=FakeTask.create_with_filemeta("hapmap"),
        sumstats_ref_task=FakeTask.create_with_filemeta("sumstats_ref"),
    )
    assert task.asset_id == "cf_dir_py"
