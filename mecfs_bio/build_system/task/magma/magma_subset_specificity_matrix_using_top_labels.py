from pathlib import Path

import pandas as pd
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.processed_gwas_data_directory_meta import (
    ProcessedGwasDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
    DataFrameWhiteSpaceSepTextFormat,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import (
    scan_dataframe,
    scan_dataframe_asset,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.magma.magma_gene_set_analysis_task import (
    GENE_SET_ANALYSIS_OUTPUT_STEM_NAME,
    MagmaGeneSetAnalysisTask,
)
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class MagmaSubsetSpecificityMatrixWithTopLabels(Task):
    """
    Task to prune a gene specificity matrix to include only the gene-property/gene-covariates that
    are found to be significant

    Intention is to mimic the operation described here:
    https://github.com/Integrative-Mental-Health-Lab/linking_cell_types_to_brain_phenotypes/blob/675b5c9b58b8762934183a3ca61ae49ad587934a/MAGMA/3.create_top_results_matrix.md

    Generally, this task is used as a preprocessing step for conditional MAGMA analysis.

    """

    _meta: Meta
    specificity_matrix_task: Task
    magma_gene_covar_analysis_task: MagmaGeneSetAnalysisTask

    nominal_sig_level: float = 0.05

    @property
    def gene_covar_result_id(self) -> AssetId:
        return self.magma_gene_covar_analysis_task.asset_id

    @property
    def spec_matrix_id(self) -> AssetId:
        return self.specificity_matrix_task.asset_id

    @property
    def spec_matrix_meta(self) -> Meta:
        return self.specificity_matrix_task.meta

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.specificity_matrix_task, self.magma_gene_covar_analysis_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        spec_matrix_asset = fetch(self.spec_matrix_id)
        spec_matrix = (
            scan_dataframe_asset(spec_matrix_asset, meta=self.spec_matrix_meta)
            .collect()
            .to_pandas()
        )

        gene_set_result_asset = fetch(self.gene_covar_result_id)
        assert isinstance(gene_set_result_asset, DirectoryAsset)

        path_to_gene_covar_result = (gene_set_result_asset.path) / (
            GENE_SET_ANALYSIS_OUTPUT_STEM_NAME + ".gsa.out"
        )

        covar_result = (
            scan_dataframe(
                path_to_gene_covar_result,
                spec=DataFrameReadSpec(
                    DataFrameWhiteSpaceSepTextFormat(comment_code="#")
                ),
            )
            .collect()
            .to_pandas()
        )
        spec_matrix_filtered = get_spec_matrix_filtered(
            covar_result=covar_result,
            spec_matrix=spec_matrix,
            nominal_sig_level=self.nominal_sig_level,
        )
        out_path = scratch_dir / self.asset_id
        spec_matrix_filtered.to_csv(out_path, index=False, sep="\t")
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        specificity_matrix_task: Task,
        magma_gene_covar_analysis_task: MagmaGeneSetAnalysisTask,
        nominal_sig_level: float = 0.05,
    ):
        source_meta = magma_gene_covar_analysis_task.meta
        assert isinstance(source_meta, ProcessedGwasDataDirectoryMeta)
        meta = ResultTableMeta(
            asset_id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            extension=".txt",
            read_spec=DataFrameReadSpec(DataFrameTextFormat(separator="\t")),
        )
        return cls(
            meta=meta,
            specificity_matrix_task=specificity_matrix_task,
            magma_gene_covar_analysis_task=magma_gene_covar_analysis_task,
            nominal_sig_level=nominal_sig_level,
        )


def get_spec_matrix_filtered(
    covar_result: pd.DataFrame,
    spec_matrix: pd.DataFrame,
    nominal_sig_level: float,
) -> pd.DataFrame:
    alpha = nominal_sig_level / len(covar_result)
    filtered_covar_result = covar_result.loc[covar_result["P"] < alpha]
    filtered_covar_result = filtered_covar_result.sort_values(by=["P"])
    sig_list = filtered_covar_result["VARIABLE"].tolist()
    cols = ["GENE"] + sig_list
    spec_matrix_filtered: pd.DataFrame = spec_matrix.loc[:, cols]
    return spec_matrix_filtered
