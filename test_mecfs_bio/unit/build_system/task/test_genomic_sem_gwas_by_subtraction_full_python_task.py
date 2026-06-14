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
    sumstats_trait,
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


def test_rejects_duplicate_aliases():
    """The two traits must have distinct aliases (they key the output columns)."""
    composite = _make_gwas_source("trait_a_data", "same")
    reference = _make_gwas_source("trait_b_data", "same")
    with pytest.raises(AssertionError, match="aliases must differ"):
        GenomicSEMGWASBySubtractionFullPythonTask.create(
            asset_id="subtraction_full_py",
            composite_trait_source=composite,
            reference_trait_source=reference,
            ld_ref_task=FakeTask.create_with_filemeta("ld_ref"),
            hapmap_snps_task=FakeTask.create_with_filemeta("hapmap"),
            sumstats_ref_task=FakeTask.create_with_filemeta("sumstats_ref"),
        )


@pytest.mark.parametrize("method", [OLS, LOGISTIC, LINEAR_PROB])
def test_sumstats_trait_carries_method(method):
    """sumstats_trait passes the source's GWAS method through unchanged."""
    src = _make_gwas_source("t_data", "t", gwas_method=method)
    df = pl.DataFrame({"SNP": ["rs1"], "P": [0.5], "effect": [0.1]})
    trait = sumstats_trait(src, df, n=10000.0)

    assert trait.gwas_method == method
    assert trait.name == "t"
    assert trait.n == 10000.0
