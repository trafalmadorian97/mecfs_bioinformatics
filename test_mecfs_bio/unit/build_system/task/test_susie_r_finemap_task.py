import tempfile
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd
import polars as pl
import polars.testing
import pytest
import rpy2.robjects as ro
import scipy.sparse
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import (
    importr,
)

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
    DataFrameWhiteSpaceSepTextFormat,
)
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.rebuilder.metadata_to_path.simple_meta_to_path import (
    SimpleMetaToPath,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.simple_hasher import (
    SimpleHasher,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.verifying_trace_info import (
    VerifyingTraceInfo,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.verifying_trace_rebuilder_core import (
    VerifyingTraceRebuilder,
)
from mecfs_bio.build_system.scheduler.topological_scheduler import topological
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.external_file_copy_task import ExternalFileCopyTask
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    ADJUSTMENT_VALUE_FILENAME,
    CS_DATA_SUBDIR,
    PIP_COLUMN,
    PIP_FILENAME,
    BroadInstituteFormatLDMatrix,
    SusieRFinemapTask,
    align_gwas_and_ld,
    extract_cs_data_tables,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
    GENE_INFO_STRAND_COL,
    BinOptions,
    HeatmapOptions,
    RegionSelectDefault,
    SusieStackPlotTask,
)
from mecfs_bio.build_system.tasks.simple_tasks import find_tasks
from mecfs_bio.build_system.wf.base_wf import SimpleWF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)


def test_align():
    """
    Verify we can correctly align GWAS data with an LD reference
    """
    gwas = pl.DataFrame(
        [
            {
                GWASLAB_CHROM_COL: 1,
                GWASLAB_POS_COL: 1,
                GWASLAB_EFFECT_ALLELE_COL: "A",
                GWASLAB_NON_EFFECT_ALLELE_COL: "C",
            },
            {
                GWASLAB_CHROM_COL: 1,
                GWASLAB_POS_COL: 2,
                GWASLAB_EFFECT_ALLELE_COL: "A",
                GWASLAB_NON_EFFECT_ALLELE_COL: "C",
            },
            {
                GWASLAB_CHROM_COL: 1,
                GWASLAB_POS_COL: 3,
                GWASLAB_EFFECT_ALLELE_COL: "A",
                GWASLAB_NON_EFFECT_ALLELE_COL: "G",
            },
        ]
    )
    reference = pl.DataFrame(
        [
            {
                GWASLAB_CHROM_COL: 1,
                GWASLAB_POS_COL: 2,
                GWASLAB_EFFECT_ALLELE_COL: "A",
                GWASLAB_NON_EFFECT_ALLELE_COL: "C",
            },
            {
                GWASLAB_CHROM_COL: 1,
                GWASLAB_POS_COL: 3,
                GWASLAB_EFFECT_ALLELE_COL: "A",
                GWASLAB_NON_EFFECT_ALLELE_COL: "G",
            },
            {
                GWASLAB_CHROM_COL: 1,
                GWASLAB_POS_COL: 4,
                GWASLAB_EFFECT_ALLELE_COL: "G",
                GWASLAB_NON_EFFECT_ALLELE_COL: "T",
            },
        ]
    )
    ld_matrix = scipy.sparse.coo_matrix(
        np.array([[1, 0.2, 0.3], [0.2, 1, 0.4], [0.3, 0.4, 1]])
    )
    rg, rr, rmat = align_gwas_and_ld(
        gwas=gwas,
        ld_labels=reference,
        ld_matrix_sparse=ld_matrix,
    )
    pl.testing.assert_frame_equal(
        rg,
        gwas[1:],
    )
    pl.testing.assert_frame_equal(rr, reference[:-1])
    np.testing.assert_array_equal(rmat, ld_matrix.toarray()[:-1, :-1])


_susie_n = 2500


@pytest.fixture(params=[[0, -1], []])
def susie_prerequisite_file_tasks(
    tmp_path: Path, request
) -> Iterator[tuple[Task, Task, Task, list[int]]]:
    n = _susie_n
    m = 100
    susie_package = importr("susieR")
    generator = np.random.default_rng(40)
    true_effects = np.zeros(m)
    causal_variants = request.param
    for cv in causal_variants:
        true_effects[cv] = 1
    q = np.linspace(1, 0, num=m)
    lamb = 0.2
    covar = (1 - lamb) * q.reshape(-1, 1) * q.reshape(1, -1) + lamb * np.eye(m)

    genotypes = generator.multivariate_normal(
        np.zeros(m), cov=covar, size=n
    )  # n rows, m columns
    phenotypes = genotypes @ true_effects.reshape(-1) + generator.normal(
        loc=0, scale=0.1, size=n
    )

    conv = ro.default_converter + pandas2ri.converter + numpy2ri.converter
    gwas_data = pd.DataFrame(
        {
            GWASLAB_CHROM_COL: [1] * m,
            GWASLAB_POS_COL: list(range(m)),
            GWASLAB_EFFECT_ALLELE_COL: "A",
            GWASLAB_NON_EFFECT_ALLELE_COL: "C",
        }
    )

    ref_data = gwas_data.copy().iloc[::-1]
    ld_matrix = np.corrcoef(genotypes.transpose())[
        ::-1, ::-1
    ]  # Flip order to ensure SUSIE task can resolve this issue

    partial_ld = ld_matrix / 2
    partial_ld_sparse = scipy.sparse.coo_matrix(partial_ld)
    with localconverter(conv):
        genotypes_r = ro.conversion.get_conversion().py2rpy(genotypes)
        phenotypes_r = ro.conversion.get_conversion().py2rpy(phenotypes)
    result = susie_package.univariate_regression(genotypes_r, phenotypes_r)
    with localconverter(conv):
        beta_hat = result.rx2("betahat")
        se_beta_hat = result.rx2("sebetahat")
    gwas_data[GWASLAB_BETA_COL] = beta_hat
    gwas_data[GWASLAB_SE_COL] = se_beta_hat
    scratch = tmp_path / "scratch"
    scratch.mkdir(parents=True, exist_ok=True)
    gwas_data_path = scratch / "gwas_data"
    ld_labels_path = scratch / "ld_labels"
    ld_matrix_path = scratch / "ld_matrix.npz"
    gwas_data.to_parquet(gwas_data_path)
    ref_data.to_parquet(ld_labels_path)
    scipy.sparse.save_npz(ld_matrix_path, partial_ld_sparse)
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        gwas_data_task = ExternalFileCopyTask(
            SimpleFileMeta(
                AssetId("gwas_data"),
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            ),
            external_path=gwas_data_path,
        )
        ld_labels_task = ExternalFileCopyTask(
            SimpleFileMeta(
                "ld_labels", read_spec=DataFrameReadSpec(DataFrameParquetFormat())
            ),
            external_path=ld_labels_path,
        )
        ld_matrix_task = ExternalFileCopyTask(
            SimpleFileMeta("ld_matrix"), external_path=ld_matrix_path
        )
        yield gwas_data_task, ld_labels_task, ld_matrix_task, causal_variants


@pytest.fixture
def dummy_ensmbl_data_task(tmp_path: Path) -> Iterator[Task]:
    dummy_data = """ENSG00000237683 1       25  50  -       AL627309.1
                    ENSG00000235249 1       77  90  +       OR4F29"""
    dummy_data_path = tmp_path / "dummy_data.txt"
    dummy_data_path.write_text(dummy_data)
    yield ExternalFileCopyTask(
        SimpleFileMeta(
            "dummy_ensmbl_data",
            read_spec=DataFrameReadSpec(
                DataFrameWhiteSpaceSepTextFormat(
                    comment_code="#",
                    col_names=[
                        "ensembl_name",
                        GENE_INFO_CHROM_COL,
                        GENE_INFO_START_COL,
                        GENE_INFO_END_COL,
                        GENE_INFO_STRAND_COL,
                        GENE_INFO_NAME_COL,
                    ],
                )
            ),
        ),
        external_path=dummy_data_path,
    )


def test_fine_mapping(
    tmp_path: Path,
    susie_prerequisite_file_tasks: tuple[Task, Task, Task, list[int]],
    dummy_ensmbl_data_task: Task,
):
    """
    Test that we can find the causal SNPs in a simple synthetic example
    """
    gwas_data_task, ld_labels_task, ld_matrix_task, causal_variants = (
        susie_prerequisite_file_tasks
    )
    susie_tsk = SusieRFinemapTask(
        meta=SimpleDirectoryMeta(AssetId("directory")),
        gwas_data_task=gwas_data_task,
        ld_labels_task=ld_labels_task,
        ld_matrix_source=BroadInstituteFormatLDMatrix(ld_matrix_task),
        effective_sample_size=_susie_n,
        max_credible_sets=10,
    )

    stack_plot_task = SusieStackPlotTask(
        meta=SimpleDirectoryMeta("susie_stack_plot"),
        susie_task=susie_tsk,
        gene_info_task=dummy_ensmbl_data_task,
        region_mode=RegionSelectDefault(),
        gene_info_pipe=IdentityPipe(),
        heatmap_options=HeatmapOptions(BinOptions(num_bins=50), mode="ld2"),
    )
    tasks = find_tasks([susie_tsk, stack_plot_task])
    wf = SimpleWF()
    info: VerifyingTraceInfo = VerifyingTraceInfo.empty()

    asset_dir = tmp_path / "asset_dir"
    asset_dir.mkdir(exist_ok=True, parents=True)
    meta_to_path = SimpleMetaToPath(root=asset_dir)

    tracer = SimpleHasher.md5_hasher()
    rebuilder = VerifyingTraceRebuilder(tracer)

    targets = [susie_tsk.asset_id, stack_plot_task.asset_id]

    # Verify that all files are created in the correct location
    store, info = topological(
        rebuilder=rebuilder,
        tasks=tasks,
        targets=targets,
        wf=wf,
        info=info,
        meta_to_path=meta_to_path,
    )
    asset = store[susie_tsk.asset_id]
    assert isinstance(asset, DirectoryAsset)
    suise_out_path = asset.path
    adjustment = pd.read_parquet(suise_out_path / ADJUSTMENT_VALUE_FILENAME)
    assert float(adjustment.iloc[0].item()) <= 0.01
    pip = pd.read_parquet(suise_out_path / PIP_FILENAME)
    for cv in causal_variants:
        assert pip[PIP_COLUMN].iloc[cv] >= 0.95
    cs_subdir_contents = list((suise_out_path / CS_DATA_SUBDIR).glob("*"))
    assert len(cs_subdir_contents) == len(causal_variants)  # 2 credible sets


def test_extract_cs_data_tables():
    """
    Test that we can correctly extract SUSIE output

    The main failure case here is when SUSIE uses its purity filters to discard the first credible set, and so
    the "SETS" object it returns does not contain L1

    This test verifies we can handle this case
    """
    gt = pl.DataFrame(
        {
            "CHROM": [20] * 5,
            "POS": [100, 200, 300, 400, 500],
            "EA": ["A", "C", "G", "T", "A"],
            "NEA": ["T", "G", "C", "A", "C"],
        }
    )

    mock_alpha = np.array(
        [
            [0.2, 0.2, 0.2, 0.2, 0.2],
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5, 0.5],
        ]
    )

    mock_mu = np.array(
        [
            [0.1, 0.1, 0.1, 0.1, 0.1],
            [0.5, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.3, 0.3],
        ]
    )

    mock_pip = np.array([1.0, 0.2, 0.2, 0.7, 0.7])

    mock_sets = {
        "cs": {
            "L2": [1],  # Refers to variant index 0
            "L3": [4, 5],  # Refers to variant indices 3 and 4
        },
        "purity": pd.DataFrame({"min.abs.corr": [1.0, 0.8]}),
    }
    cs_table_dict = extract_cs_data_tables(
        alpha=mock_alpha,
        mu=mock_mu,
        pip=mock_pip,
        sets=mock_sets,
        gt=gt,
    )
    combined_cs = pl.concat(cs_table_dict.values(), how="vertical")
    l2_mu = combined_cs.filter(pl.col("cs") == "L2")["mu"][0]
    assert l2_mu == 0.5
