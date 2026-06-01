"""
Unit tests for GenomicSEMGWASBySubtractionFullPythonTask -- the fully-Python
(no rpy2) GWAS-by-subtraction task. Tests cover Python-side wiring only (deps,
meta, k=2 assertion, the gwas_method -> SumstatsTrait flag mapping). The
numerical correctness of each stage lives in test_genomic_sem_ldsc.py,
test_genomic_sem_munge.py, test_genomic_sem_sumstats.py, and
test_gwas_by_subtraction_kernel.py.
"""

import polars as pl
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
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_config import (
    LINEAR_PROB,
    LOGISTIC,
    OLS,
    GenomicSEMGWASSumstatsSource,
    GenomicSEMSumstatsSource,
    GWASMethod,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_gwas_by_subtraction_full_python_task import (
    GenomicSEMGWASBySubtractionFullPythonTask,
    _sumstats_trait,
)


def _make_gwas_source(
    asset_id: str, alias: str, gwas_method: GWASMethod = OLS
) -> GenomicSEMGWASSumstatsSource:
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
    return GenomicSEMGWASSumstatsSource(source=inner, gwas_method=gwas_method)


def test_create_meta_and_deps():
    ld_ref = FakeTask.create_with_filemeta("ld_ref")
    hapmap = FakeTask.create_with_filemeta("hapmap")
    ref = FakeTask.create_with_filemeta("sumstats_ref")
    src_a = _make_gwas_source("trait_a_data", "a")
    src_b = _make_gwas_source("trait_b_data", "b")

    task = GenomicSEMGWASBySubtractionFullPythonTask.create(
        asset_id="subtraction_full_py",
        sources=[src_a, src_b],
        ld_ref_task=ld_ref,
        hapmap_snps_task=hapmap,
        sumstats_ref_task=ref,
    )

    assert isinstance(task.meta, ResultDirectoryMeta)
    assert task.meta.asset_id == "subtraction_full_py"
    assert task.meta.project == "genomic_sem"
    assert task.deps == [ld_ref, hapmap, ref, src_a.task, src_b.task]


def test_requires_exactly_two_sources():
    src_a = _make_gwas_source("trait_a_data", "a")
    with pytest.raises(AssertionError):
        GenomicSEMGWASBySubtractionFullPythonTask.create(
            asset_id="subtraction_full_py",
            sources=[src_a],
            ld_ref_task=FakeTask.create_with_filemeta("ld_ref"),
            hapmap_snps_task=FakeTask.create_with_filemeta("hapmap"),
            sumstats_ref_task=FakeTask.create_with_filemeta("sumstats_ref"),
        )


def test_rejects_more_than_two_traits():
    src_a = _make_gwas_source("trait_a_data", "a")
    src_b = _make_gwas_source("trait_b_data", "b")
    src_c = _make_gwas_source("trait_c_data", "c")
    with pytest.raises(AssertionError, match="exactly 2 traits"):
        GenomicSEMGWASBySubtractionFullPythonTask.create(
            asset_id="subtraction_full_py",
            sources=[src_a, src_b, src_c],
            ld_ref_task=FakeTask.create_with_filemeta("ld_ref"),
            hapmap_snps_task=FakeTask.create_with_filemeta("hapmap"),
            sumstats_ref_task=FakeTask.create_with_filemeta("sumstats_ref"),
        )


def test_uses_simple_directory_meta_when_passed():
    src_a = _make_gwas_source("a_data", "a")
    src_b = _make_gwas_source("b_data", "b")
    task = GenomicSEMGWASBySubtractionFullPythonTask(
        meta=SimpleDirectoryMeta(AssetId("subtraction_dir_py")),
        sources=[src_a, src_b],
        ld_ref_task=FakeTask.create_with_filemeta("ld_ref"),
        hapmap_snps_task=FakeTask.create_with_filemeta("hapmap"),
        sumstats_ref_task=FakeTask.create_with_filemeta("sumstats_ref"),
    )
    assert task.asset_id == "subtraction_dir_py"


@pytest.mark.parametrize(
    "method, expected",
    [
        (OLS, (True, False, False)),
        (LOGISTIC, (False, True, False)),
        (LINEAR_PROB, (False, False, True)),
    ],
)
def test_sumstats_trait_flag_mapping(method, expected):
    """gwas_method maps onto exactly one of (ols, se_logit, linprob)."""
    src = _make_gwas_source("t_data", "t", gwas_method=method)
    df = pl.DataFrame({"SNP": ["rs1"], "P": [0.5], "effect": [0.1]})
    trait = _sumstats_trait(src, df, n=10000.0)

    assert (trait.ols, trait.se_logit, trait.linprob) == expected
    assert trait.name == "t"
    assert trait.n == 10000.0
