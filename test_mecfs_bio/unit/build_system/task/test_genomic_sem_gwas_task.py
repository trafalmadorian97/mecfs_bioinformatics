"""
Unit tests for GenomicSEMCommonFactorGWASTask and GenomicSEMUserGWASTask.

Covers Python-side wiring: dependency graph, source/alias validation,
sub-component handling, GWAS-method-to-flags conversion, filename
sanitisation, and config defaults. The actual rpy2/R workflow against a real
LD reference and 1000G panel is left to a future system test.
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
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_gwas_task import (
    COMMON_FACTOR_GWAS_FILENAME,
    GWAS_RESULTS_SUBDIR,
    LINEAR_PROB,
    LOGISTIC,
    OLS,
    GenomicSEMCommonFactorGWASTask,
    GenomicSEMGWASRunConfig,
    GenomicSEMGWASSumstatsSource,
    GenomicSEMSumstatsConfig,
    GenomicSEMUserGWASTask,
    GWASMethod,
    _gwas_method_flags,
    _sanitize_component_name,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_task import (
    GenomicSEMSumstatsSource,
)


def _make_gwas_source(
    asset_id: str,
    alias: str,
    sample_info,
    gwas_method: GWASMethod = OLS,
) -> GenomicSEMGWASSumstatsSource:
    inner_task = FakeTask(
        SimpleFileMeta(
            AssetId(asset_id),
            read_spec=DataFrameReadSpec(DataFrameTextFormat(",")),
        ),
    )
    inner = GenomicSEMSumstatsSource(
        task=inner_task, alias=alias, sample_info=sample_info
    )
    return GenomicSEMGWASSumstatsSource(source=inner, gwas_method=gwas_method)


# ---- GenomicSEMCommonFactorGWASTask ----------------------------------------


def test_common_factor_create_meta_and_deps():
    """
    Factory should construct ResultDirectoryMeta and deps should include
    LD ref, hapmap, sumstats ref, and each source task.
    """
    ld_ref = FakeTask.create_with_filemeta("ld_ref")
    hapmap = FakeTask.create_with_filemeta("hapmap")
    ref = FakeTask.create_with_filemeta("sumstats_ref")
    src_a = _make_gwas_source(
        "trait_a_data", "a", QuantPhenotype(total_sample_size=10000)
    )
    src_b = _make_gwas_source(
        "trait_b_data", "b", QuantPhenotype(total_sample_size=20000)
    )

    task = GenomicSEMCommonFactorGWASTask.create(
        asset_id="cf_gwas",
        sources=[src_a, src_b],
        ld_ref_task=ld_ref,
        hapmap_snps_task=hapmap,
        sumstats_ref_task=ref,
    )

    assert isinstance(task.meta, ResultDirectoryMeta)
    assert task.meta.asset_id == "cf_gwas"
    assert task.meta.project == "genomic_sem"
    assert task.deps == [ld_ref, hapmap, ref, src_a.task, src_b.task]


def test_common_factor_requires_at_least_two_sources():
    src_a = _make_gwas_source(
        "trait_a_data", "a", QuantPhenotype(total_sample_size=10000)
    )
    with pytest.raises(AssertionError):
        GenomicSEMCommonFactorGWASTask.create(
            asset_id="cf_gwas",
            sources=[src_a],
            ld_ref_task=FakeTask.create_with_filemeta("ld_ref"),
            hapmap_snps_task=FakeTask.create_with_filemeta("hapmap"),
            sumstats_ref_task=FakeTask.create_with_filemeta("sumstats_ref"),
        )


def test_common_factor_rejects_duplicate_aliases():
    src_a = _make_gwas_source(
        "trait_a_data", "a", QuantPhenotype(total_sample_size=10000)
    )
    src_b = _make_gwas_source(
        "trait_b_data", "a", QuantPhenotype(total_sample_size=20000)
    )
    with pytest.raises(AssertionError):
        GenomicSEMCommonFactorGWASTask.create(
            asset_id="cf_gwas",
            sources=[src_a, src_b],
            ld_ref_task=FakeTask.create_with_filemeta("ld_ref"),
            hapmap_snps_task=FakeTask.create_with_filemeta("hapmap"),
            sumstats_ref_task=FakeTask.create_with_filemeta("sumstats_ref"),
        )


def test_common_factor_uses_simple_directory_meta_when_passed():
    """Direct construction with any Meta should be accepted (used by tests)."""
    src_a = _make_gwas_source("a_data", "a", QuantPhenotype(total_sample_size=10000))
    src_b = _make_gwas_source("b_data", "b", QuantPhenotype(total_sample_size=10000))
    task = GenomicSEMCommonFactorGWASTask(
        meta=SimpleDirectoryMeta(AssetId("cf_dir")),
        sources=[src_a, src_b],
        ld_ref_task=FakeTask.create_with_filemeta("ld_ref"),
        hapmap_snps_task=FakeTask.create_with_filemeta("hapmap"),
        sumstats_ref_task=FakeTask.create_with_filemeta("sumstats_ref"),
    )
    assert task.asset_id == "cf_dir"


# ---- GenomicSEMUserGWASTask ------------------------------------------------


def test_user_gwas_create_meta_and_deps():
    ld_ref = FakeTask.create_with_filemeta("ld_ref")
    hapmap = FakeTask.create_with_filemeta("hapmap")
    ref = FakeTask.create_with_filemeta("sumstats_ref")
    src_a = _make_gwas_source(
        "trait_a_data", "a", QuantPhenotype(total_sample_size=10000)
    )
    src_b = _make_gwas_source(
        "trait_b_data", "b", QuantPhenotype(total_sample_size=20000)
    )

    task = GenomicSEMUserGWASTask.create(
        asset_id="user_gwas",
        sources=[src_a, src_b],
        ld_ref_task=ld_ref,
        hapmap_snps_task=hapmap,
        sumstats_ref_task=ref,
        factor_model="F1 =~ a + b\nF1 ~ SNP",
        sub_components=["F1~SNP"],
    )

    assert isinstance(task.meta, ResultDirectoryMeta)
    assert task.deps == [ld_ref, hapmap, ref, src_a.task, src_b.task]
    assert task.factor_model == "F1 =~ a + b\nF1 ~ SNP"
    assert list(task.sub_components) == ["F1~SNP"]


def test_user_gwas_requires_at_least_two_sources():
    src_a = _make_gwas_source(
        "trait_a_data", "a", QuantPhenotype(total_sample_size=10000)
    )
    with pytest.raises(AssertionError):
        GenomicSEMUserGWASTask.create(
            asset_id="user_gwas",
            sources=[src_a],
            ld_ref_task=FakeTask.create_with_filemeta("ld_ref"),
            hapmap_snps_task=FakeTask.create_with_filemeta("hapmap"),
            sumstats_ref_task=FakeTask.create_with_filemeta("sumstats_ref"),
            factor_model="F1 =~ a",
            sub_components=["F1~SNP"],
        )


def test_user_gwas_requires_at_least_one_sub_component():
    src_a = _make_gwas_source(
        "trait_a_data", "a", QuantPhenotype(total_sample_size=10000)
    )
    src_b = _make_gwas_source(
        "trait_b_data", "b", QuantPhenotype(total_sample_size=20000)
    )
    with pytest.raises(AssertionError):
        GenomicSEMUserGWASTask.create(
            asset_id="user_gwas",
            sources=[src_a, src_b],
            ld_ref_task=FakeTask.create_with_filemeta("ld_ref"),
            hapmap_snps_task=FakeTask.create_with_filemeta("hapmap"),
            sumstats_ref_task=FakeTask.create_with_filemeta("sumstats_ref"),
            factor_model="F1 =~ a + b\nF1 ~ SNP",
            sub_components=[],
        )


def test_user_gwas_rejects_colliding_sub_component_filenames():
    """
    Two sub_components that sanitise to the same filename are a configuration
    error — the second would silently overwrite the first.
    """
    src_a = _make_gwas_source(
        "trait_a_data", "a", QuantPhenotype(total_sample_size=10000)
    )
    src_b = _make_gwas_source(
        "trait_b_data", "b", QuantPhenotype(total_sample_size=20000)
    )
    with pytest.raises(AssertionError):
        GenomicSEMUserGWASTask.create(
            asset_id="user_gwas",
            sources=[src_a, src_b],
            ld_ref_task=FakeTask.create_with_filemeta("ld_ref"),
            hapmap_snps_task=FakeTask.create_with_filemeta("hapmap"),
            sumstats_ref_task=FakeTask.create_with_filemeta("sumstats_ref"),
            factor_model="F1 =~ a + b\nF1 ~ SNP",
            sub_components=["F1~SNP", "F1=~SNP"],  # both → 'F1_SNP'
        )


# ---- helpers ---------------------------------------------------------------


def test_gwas_method_flags_one_per_trait():
    """
    Each source's GWASMethod produces exactly one TRUE in the
    (se_logit, OLS, linprob) triple, in the slot matching its method.
    """
    src_ols = _make_gwas_source("a", "a", QuantPhenotype(total_sample_size=1), OLS)
    src_log = _make_gwas_source("b", "b", QuantPhenotype(total_sample_size=1), LOGISTIC)
    src_lin = _make_gwas_source(
        "c", "c", QuantPhenotype(total_sample_size=1), LINEAR_PROB
    )

    se_logit, ols, linprob = _gwas_method_flags([src_ols, src_log, src_lin])

    assert se_logit == [False, True, False]
    assert ols == [True, False, False]
    assert linprob == [False, False, True]


def test_sanitize_component_name():
    assert _sanitize_component_name("F1~SNP") == "F1_SNP"
    assert _sanitize_component_name("F1=~a") == "F1_a"
    assert _sanitize_component_name("a~~b") == "a_b"
    assert _sanitize_component_name("~F1") == "F1"


# ---- config defaults -------------------------------------------------------


def test_sumstats_config_defaults():
    """sumstats() info filter default is 0.6 in GenomicSEM, lower than munge's 0.9."""
    cfg = GenomicSEMSumstatsConfig()
    assert cfg.info_filter == pytest.approx(0.6)
    assert cfg.maf_filter == pytest.approx(0.01)
    assert cfg.keep_indel is False
    assert cfg.parallel is False
    assert cfg.cores is None
    assert cfg.ambig is False


def test_run_config_defaults():
    cfg = GenomicSEMGWASRunConfig()
    assert cfg.estimation == "DWLS"
    assert cfg.parallel is True
    assert cfg.cores is None
    assert cfg.gc_correction == "standard"
    assert cfg.toler is False
    assert cfg.snp_se is False
    assert cfg.smooth_check is False


def test_output_filename_constants():
    """Sanity check on the constants downstream code may key on."""
    assert GWAS_RESULTS_SUBDIR == "gwas_results"
    assert COMMON_FACTOR_GWAS_FILENAME.endswith(".parquet")
