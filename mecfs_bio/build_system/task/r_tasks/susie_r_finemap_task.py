"""
Finemap a GWAS locus using SUSIE
"""

from rich.pretty import pprint
from pathlib import Path, PurePath

import numpy as np
import pandas as pd
import polars as pl
import rpy2.robjects as ro
import scipy.sparse
import structlog
import yaml
from attrs import frozen
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import (
    importr,
)
from scipy.sparse import coo_matrix, csr_matrix

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.two_sample_mr_task import RPackageType
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)

logger = structlog.get_logger()


@frozen
class BroadInstituteFormatLDMatrix:
    """
    Wraps a task to recovery an LD matrix stored in Broad Institute Format.
    This means that only the lower triangular portion is stored, and the diagonal is halved.
    To recover the full LD matrix, compute M +M^T
    """

    task: Task


LDMatrixSource = BroadInstituteFormatLDMatrix

KRIGING_PLOT_FILENAME = "kriging_diagnostic_plot.png"
ADJUSTMENT_VALUE_FILENAME = "adjustment_value.parquet"

ALPHA_FILENAME = "alpha.parquet"
PIP_FILENAME = "pip.parquet"
MU_FILENAME = "mu.parquet"
MU2_FILENAME = "mu2.parquet"
PURITY_FILENAME = "purity.parquet"
SETS_FILENAME = "sets.yaml"
COMBINED_CS_FILENAME = "combined_cs.parquet"
FILTERED_GWAS_FILENAME = "filtered_gwas.parquet"
FILTERED_LD_FILENAME = "filtered_ld.npy"

CS_DATA_SUBDIR = "credible_set_data"

PIP_COLUMN = "PIP"
ALPHA_COLUMN_NAME = "alpha"
MU_COLUMN_NAME = "mu"
MU2_COLUMN_NAME = "mu2"
CS_COLUMN = "cs"


@frozen
class SusieRFinemapTask(Task):
    """
    This task uses susie_rss function from the R library susieR to fine map a GWAS at a locus.
    It requires an LD matrix from a closely matched reference panel
    """

    _meta: Meta
    gwas_data_task: Task
    ld_labels_task: Task
    ld_matrix_source: LDMatrixSource
    effective_sample_size: float
    gwas_data_pipe: DataProcessingPipe = IdentityPipe()
    ld_labels_pipe: DataProcessingPipe = IdentityPipe()
    subsample: int | None = None
    max_credible_sets: int = 10

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.gwas_data_task, self.ld_labels_task, self.ld_matrix_source.task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        susie_package = importr("susieR")
        ggplot2 = importr("ggplot2")
        conv = ro.default_converter + pandas2ri.converter + numpy2ri.converter

        gwas_asset = fetch(self.gwas_data_task.asset_id)
        gwas_table = (
            self.gwas_data_pipe.process(
                scan_dataframe_asset(gwas_asset, self.gwas_data_task.meta)
            )
            .collect()
            .to_polars()
        )

        ld_labels_asset = fetch(self.ld_labels_task.asset_id)
        ld_labels_table = (
            self.ld_labels_pipe.process(
                scan_dataframe_asset(ld_labels_asset, self.ld_labels_task.meta)
            )
            .collect()
            .to_polars()
        )

        ld_matrix_asset = fetch(self.ld_matrix_source.task.asset_id)
        assert isinstance(ld_matrix_asset, FileAsset)
        ld_matrix_sparse = _load_ld_matrix(
            path=ld_matrix_asset.path, ld_labels_table=ld_labels_table
        )

        gwas_table, ld_labels_table, ld_matrix = align_gwas_and_ld(
            gwas=gwas_table,
            ld_labels=ld_labels_table,
            ld_matrix_sparse=ld_matrix_sparse,
        )
        del ld_matrix_sparse
        gwas_table, ld_labels_table, ld_matrix = apply_subsample(
            gwas_table, ld_labels_table, ld_matrix, subsample=self.subsample
        )
        assert len(gwas_table) == len(ld_labels_table) == len(ld_matrix)
        logger.debug(f"Dimensions of LD matrix to use: {ld_matrix.shape}")

        ld_matrix = make_psd_corr(ld_matrix)

        zscores_pandas = (
            gwas_table[GWASLAB_BETA_COL] / gwas_table[GWASLAB_SE_COL]
        ).to_pandas()

        with localconverter(conv):
            zscores_r = ro.conversion.get_conversion().py2rpy(zscores_pandas)
            ld_matrix_r = ro.conversion.get_conversion().py2rpy(ld_matrix)
            logger.debug("Running Summary Statistics/ LD matrix diagnostics")
            adjustment = susie_package.estimate_s_rss(
                zscores_r, ld_matrix_r, n=int(self.effective_sample_size)
            )
            logger.debug(f"LD matrix adjustment parameter is: {adjustment}")
        _make_diagnostic_plot(
            zscores_r=zscores_r,
            ld_matrix_r=ld_matrix_r,
            effective_sample_size=int(self.effective_sample_size),
            susie_package=susie_package,
            ggplot2_package=ggplot2,
            scratch_dir=scratch_dir,
        )
        ld_matrix = (1 - adjustment) * ld_matrix + adjustment * np.eye(
            ld_matrix.shape[0]
        )
        check_symmetric(ld_matrix)
        with localconverter(conv):
            ld_matrix_r = ro.conversion.get_conversion().py2rpy(ld_matrix)
        _save_adjustment(adjustment=float(adjustment), scratch_dir=scratch_dir)

        logger.debug("Running SUSIE")
        susie_result = susie_package.susie_rss(
            zscores_r,
            ld_matrix_r,
            n=int(self.effective_sample_size),
            L=self.max_credible_sets,
        )
        py_result = convert_r_named_list_to_python_dict(susie_result)
        _check_converged(py_result)
        write_result(
            scratch_dir,
            py_result=py_result,
            gwas_table=gwas_table,
            ld_matrix=ld_matrix,
            L=self.max_credible_sets,
        )

        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        gwas_data_task: Task,
        ld_labels_task: Task,
        ld_matrix_source: LDMatrixSource,
        effective_sample_size: float,
        gwas_data_pipe: DataProcessingPipe = IdentityPipe(),
        ld_labels_pipe: DataProcessingPipe = IdentityPipe(),
        subsample: int | None = None,
        max_credible_sets: int = 10,
    ):
        source_meta = gwas_data_task.meta
        meta: Meta
        if isinstance(source_meta, FilteredGWASDataMeta):
            meta = ResultDirectoryMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                sub_dir=PurePath("analysis"),
            )
        else:
            raise NotImplementedError(f"No handler for meta {source_meta} implemented")
        return cls(
            meta=meta,
            gwas_data_task=gwas_data_task,
            ld_labels_task=ld_labels_task,
            ld_matrix_source=ld_matrix_source,
            effective_sample_size=effective_sample_size,
            gwas_data_pipe=gwas_data_pipe,
            ld_labels_pipe=ld_labels_pipe,
            max_credible_sets=max_credible_sets,
            subsample=subsample,
        )


def align_gwas_and_ld(
    gwas: pl.DataFrame, ld_labels: pl.DataFrame, ld_matrix_sparse: coo_matrix
) -> tuple[pl.DataFrame, pl.DataFrame, np.ndarray]:
    """
    Slice the reference LD matrix and the GWAS data so that they only include genetic variants in their intersection
    """
    gwas = gwas.sort(by=GWASLAB_POS_COL)
    gwas = gwas.with_row_index(name="gwas_index")
    ld_labels = ld_labels.with_row_index(name="ld_index")
    joined = gwas.join(
        ld_labels,
        on=[
            GWASLAB_CHROM_COL,
            GWASLAB_POS_COL,
            GWASLAB_EFFECT_ALLELE_COL,
            GWASLAB_NON_EFFECT_ALLELE_COL,
        ],maintain_order="left"
    )
    ld_matrix = csr_matrix(ld_matrix_sparse).toarray()
    ld_matrix = ld_matrix[
        joined["ld_index"].to_numpy().reshape(-1, 1),
        joined["ld_index"].to_numpy().reshape(1, -1),
    ]
    return (
        gwas[joined["gwas_index"]].drop("gwas_index"),
        ld_labels[joined["ld_index"]].drop("ld_index"),
        ld_matrix,
    )


def make_psd_corr(matrix: np.ndarray, tol: float = 1e-4) -> np.ndarray:
    """
    Make a correlation matrix positive semidefinite

    """
    logger.debug("ensuring LD matrix is positive definite")
    eigs, eigvecs = np.linalg.eigh(matrix)
    if (eigs >= tol).all():
        logger.debug(f"eigenvalues exceed {tol}.  No adjustment necessary.")
        return matrix
    logger.debug(f"Smallest eigenvalue is {np.min(eigs)} < {tol}.  Applying adjustment")

    eigs[eigs < tol] = tol
    reconstructed = eigvecs @ np.diag(eigs) @ eigvecs.T
    d = np.diag(reconstructed)
    reconstructed = reconstructed / np.sqrt(d.reshape(-1, 1) * d.reshape(1, -1))
    assert abs(np.diag(reconstructed) - 1).max() <= 1e-5
    return reconstructed


def apply_subsample(
    gwas_table: pl.DataFrame,
    ld_labels: pl.DataFrame,
    ld_matrix: np.ndarray,
    subsample: int | None,
) -> tuple[pl.DataFrame, pl.DataFrame, np.ndarray]:
    """
    Subsample gwasdata
    """
    if subsample is None:
        return gwas_table, ld_labels, ld_matrix
    logger.debug("Subsampling...")
    return (
        gwas_table[::subsample],
        ld_labels[::subsample],
        ld_matrix[::subsample, ::subsample],
    )


def convert_r_named_list_to_python_dict(rdata) -> dict:
    return dict(zip(rdata.names, list(rdata)))


def r_to_py_recursive(r_obj):
    conv = ro.default_converter + pandas2ri.converter + numpy2ri.converter
    if isinstance(r_obj, ro.DataFrame):
        with localconverter(conv):
            alpha = ro.conversion.get_conversion().rpy2py(r_obj)
        return pandas2ri.rpy2py(r_obj)

    if isinstance(r_obj, ro.ListVector):
        names = r_obj.names

        if names is None:
            return [r_to_py_recursive(item) for item in r_obj]

        return {name: r_to_py_recursive(r_obj.rx2(name)) for name in names}

    try:
        if hasattr(r_obj, "__len__"):
            if len(r_obj) == 1:
                return r_obj[0]
            return list(r_obj)
    except:
        pass

    return r_obj


def _check_converged(py_result: dict):
    conv = ro.default_converter
    with localconverter(conv):
        converged = ro.conversion.get_conversion().rpy2py(py_result["converged"])

    logger.debug(f"Convergence: {converged}")
    assert bool(converged), "SUSIE failed to converge"


def write_result(
    directory: Path, py_result: dict, gwas_table: pl.DataFrame, L: int,
        ld_matrix: np.ndarray,
) -> None:
    gt = gwas_table.select(
        [
            GWASLAB_CHROM_COL,
            GWASLAB_POS_COL,
            GWASLAB_EFFECT_ALLELE_COL,
            GWASLAB_NON_EFFECT_ALLELE_COL,
        ]
    )
    credible_set_names = [f"L{i}" for i in range(1, L + 1)]
    variant_names = [
        f"ch{chrom}_{pos}_{ea}_{nea}"
        for chrom, pos, ea, nea in zip(
            gt[GWASLAB_CHROM_COL],
            gt[GWASLAB_POS_COL],
            gt[GWASLAB_EFFECT_ALLELE_COL],
            gt[GWASLAB_NON_EFFECT_ALLELE_COL],
        )
    ]

    conv = ro.default_converter + pandas2ri.converter + numpy2ri.converter

    sets = r_to_py_recursive(py_result["sets"])
    with localconverter(conv):
        alpha = ro.conversion.get_conversion().rpy2py(py_result["alpha"])
        pip = ro.conversion.get_conversion().rpy2py(py_result["pip"])
        mu = ro.conversion.get_conversion().rpy2py(py_result["mu"])
        mu2 = ro.conversion.get_conversion().rpy2py(py_result["mu2"])
    pd.DataFrame(alpha, columns=variant_names, index=credible_set_names).to_parquet(
        directory / ALPHA_FILENAME
    )
    pd.DataFrame(pip, index=variant_names, columns=[PIP_COLUMN]).to_parquet(
        directory / PIP_FILENAME
    )
    pd.DataFrame(mu, columns=variant_names, index=credible_set_names).to_parquet(
        directory / MU_FILENAME
    )
    pd.DataFrame(mu2, columns=variant_names, index=credible_set_names).to_parquet(
        directory / MU2_FILENAME
    )
    pd.DataFrame(sets["purity"]).to_parquet(directory / PURITY_FILENAME)
    sets.pop("purity")
    with open(directory / "sets.yaml", "w") as f:
        yaml.dump(
            sets,
            f,
            sort_keys=True,
            default_flow_style=False,
            allow_unicode=True,
        )
    subdir = directory / CS_DATA_SUBDIR
    subdir.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Found {len(sets['cs'])} credible sets")
    cs_data_tables = []
    for set_number, (set_name, set_values) in enumerate(sets["cs"].items()):
        set_values = (
            np.array(set_values).reshape(-1) - 1
        )  # Since r uses 1-based indexing
        gwas_for_set: pl.DataFrame = gt[set_values]
        alpha_for_set = pl.Series(
            values=alpha[set_number, set_values].reshape(-1), name=ALPHA_COLUMN_NAME
        )
        mu_for_set = pl.Series(
            values=mu[set_number, set_values].reshape(-1), name=MU_COLUMN_NAME
        )
        # mu2_for_set = pl.Series(values=mu2[set_number,set_values].reshape(-1), name=MU2_COLUMN_NAME)
        pip_for_set = pl.Series(values=pip[set_values].reshape(-1), name=PIP_COLUMN)
        gwas_for_set = gwas_for_set.with_columns(alpha_for_set, mu_for_set, pip_for_set)
        logger.debug(
            f"Credible set {set_name}\n {gwas_for_set.sort(by=PIP_COLUMN, descending=True)}"
        )
        gwas_for_set.write_parquet(subdir / f"{set_name}.parquet")
        cs_data_tables.append(
            gwas_for_set.with_columns(pl.lit(set_name).alias(CS_COLUMN))
        )
    pl.concat(cs_data_tables,how="vertical").write_parquet(directory/COMBINED_CS_FILENAME)
    gwas_table.write_parquet(directory / FILTERED_GWAS_FILENAME)
    np.save(directory/FILTERED_LD_FILENAME, ld_matrix)




def check_symmetric(array: np.ndarray, tol=1e-6):
    assert np.max(np.abs(array - array.T)) <= tol


def _save_adjustment(adjustment: float, scratch_dir: Path):
    adjustment_df = pd.DataFrame(
        {"Adjustment": [float(adjustment)]},
    )
    adjustment_df.to_parquet(scratch_dir / ADJUSTMENT_VALUE_FILENAME)


def _load_ld_matrix(path: Path, ld_labels_table: pl.DataFrame) -> coo_matrix:
    logger.debug(f"loading ld matrix from {path}")
    partial_ld_matrix = scipy.sparse.load_npz(path)
    ld_matrix = partial_ld_matrix + partial_ld_matrix.transpose()
    logger.debug("done loading ld matrix")

    assert ld_matrix.shape[0] == len(ld_labels_table)
    assert (abs(ld_matrix.diagonal() - 1)).max() <= 1e-4
    return coo_matrix(ld_matrix)


def _make_diagnostic_plot(
    zscores_r,
    ld_matrix_r,
    effective_sample_size: int,
    susie_package: RPackageType,
    ggplot2_package: RPackageType,
    scratch_dir: Path,
):
    logger.debug("Creating diagnostic plot")
    plot_and_table = susie_package.kriging_rss(
        zscores_r, ld_matrix_r, n=effective_sample_size
    )
    r_plot_object = plot_and_table.rx2("plot")
    ggplot2_package.ggsave(
        filename=str(scratch_dir / KRIGING_PLOT_FILENAME),
        plot=r_plot_object,
        width=7,
        height=7,
        dpi=300,
        device="png",
    )
