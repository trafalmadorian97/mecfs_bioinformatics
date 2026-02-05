from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import polars.testing
import rpy2.robjects as ro
import scipy.sparse
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import (
    importr,
)

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    ADJUSTMENT_VALUE_FILENAME,
    CS_DATA_SUBDIR,
    PIP_COLUMN,
    PIP_FILENAME,
    BroadInstituteFormatLDMatrix,
    SusieRFinemapTask,
    align_gwas_and_ld,
)
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
        ld_matrix=ld_matrix,
    )
    pl.testing.assert_frame_equal(
        rg,
        gwas[1:],
    )
    pl.testing.assert_frame_equal(rr, reference[:-1])
    np.testing.assert_array_equal(rmat, ld_matrix.toarray()[:-1, :-1])


def test_fine_mapping(tmp_path: Path):
    """
    Test that we can find the causal SNPs in a simple synthetic example
    """
    n = 2500
    m = 100
    susie_package = importr("susieR")
    generator = np.random.default_rng(40)
    true_effects = np.zeros(m)
    true_effects[0] = 1
    true_effects[-1] = 1
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
    tsk = SusieRFinemapTask(
        meta=SimpleDirectoryMeta(AssetId("directory")),
        gwas_data_task=FakeTask(
            SimpleFileMeta(
                AssetId("gwas_data"),
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            ),
        ),
        ld_labels_task=FakeTask(
            SimpleFileMeta(
                "ld_labels", read_spec=DataFrameReadSpec(DataFrameParquetFormat())
            )
        ),
        ld_matrix_source=BroadInstituteFormatLDMatrix(
            FakeTask(SimpleFileMeta("ld_matrix")),
        ),
        effective_sample_size=n,
        max_credible_sets=10,
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "gwas_data":
            return FileAsset(gwas_data_path)
        if asset_id == "ld_labels":
            return FileAsset(ld_labels_path)
        if asset_id == "ld_matrix":
            return FileAsset(ld_matrix_path)
        raise ValueError("unknown id")

    tsk.execute(scratch_dir=scratch, fetch=fetch, wf=SimpleWF())
    adjustment = pd.read_parquet(scratch / ADJUSTMENT_VALUE_FILENAME)
    assert float(adjustment.iloc[0].item()) <= 0.01
    pip = pd.read_parquet(scratch / PIP_FILENAME)
    assert pip[PIP_COLUMN].iloc[0] >= 0.95
    assert pip[PIP_COLUMN].iloc[99] >= 0.95
    cs_subdir_contents = list((scratch / CS_DATA_SUBDIR).glob("*"))
    assert len(cs_subdir_contents) == 2  # 2 credible sets
