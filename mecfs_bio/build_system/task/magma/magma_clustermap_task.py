from typing import List
import structlog
logger = structlog.get_logger()

import pandas as pd
import seaborn as sns
from pathlib import Path

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.plot_meta import GWASPlotDirectoryMeta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.magma.magma_gene_analysis_task import MagmaGeneAnalysisTask, \
    read_magma_gene_analysis_result
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.magma_constants import MAGMA_GENE_COL, MAGMA_P_COLUMN, MAGMA_Z_COLUMN


# to try:
# set a minimum value for the color scale, since we are doing a one-sided test in MAGMA.
#  or Choose a divergent color scale..
# Change normalization to be consistent with regression: subtract mean without dividing by std.
# or divide columns by std
# or COMPUTE COLUMN-wise Z-score across entire matrix, not just selected genes


##HMM
# What magma says:
# Pattern of gene significance matches pattern of gene expression in brain

# What I am looking for: Significance genes are highly expressed in brain.


@frozen
class TopGenesMode:
    num_genes: int
    source_task: MagmaGeneAnalysisTask
    p_value_thresh:float

@frozen
class GeneNamesMode:
    gene_names: List[str]
    gene_thesaurus: Task


GeneSelectMode = TopGenesMode | GeneNamesMode


@frozen
class TopClustersMode:
    num_clusters: int
    source_task: Task

@frozen
class AllClustersMode:
    pass

ClusterSelectMode = TopClustersMode | AllClustersMode


@frozen
class AuxInfoFrameSpec:
    aux_info_task:Task
    col_in_aux_info:str


@frozen
class GeneWiseSpecificityNormalization:
    pass


@frozen
class ClusterWiseZScoreNormalization:
    pass


NormalizationStrategy = GeneWiseSpecificityNormalization | ClusterWiseZScoreNormalization




@frozen
class VlimForGeneWiseSpec:
    max_multiple:float

Vlimspec = VlimForGeneWiseSpec




@frozen
class ExpressionMatrixClusterMapTask(Task):
    """
    Task to make a matrix plot with genes as rows and tissue/cell types as columns by subsetting a MAGMA source matrix

    Various possible gene subsetting rules:
    - Most significant genes
    - Most significant "independent" genes?

    Various possible tissue/cell type rules
    - Most significant tissue/cell types
    - Most significant "independent" tissue/cell types
    """
    _meta: Meta
    specificity_matrix_task: Task
    gene_select_mode:GeneSelectMode
    cluster_select_mode: ClusterSelectMode
    gene_col_name: str
    normalizations: list[NormalizationStrategy]
    plot_params: dict
    vlimspec: Vlimspec|None
    aux_info_spec: AuxInfoFrameSpec | None = None



    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        result = [self.specificity_matrix_task]
        if isinstance(self.gene_select_mode, TopGenesMode):
            result.append(self.gene_select_mode.source_task)
        if isinstance(self.cluster_select_mode, TopClustersMode):
            result.append(self.cluster_select_mode.source_task)
        if self.aux_info_spec is not None:
            result.append(self.aux_info_spec.aux_info_task)
        return result

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        gene_cluster_raw_asset = fetch(self.specificity_matrix_task.asset_id)

        gene_cluster_matrix_raw = scan_dataframe_asset(
            gene_cluster_raw_asset,
            meta=self.specificity_matrix_task.meta,
        ).collect().to_pandas()
        _debugging(
            gene_select_mode=self.gene_select_mode,
            fetch=fetch,
            gene_cluster_matrix_raw=gene_cluster_matrix_raw,
        )

        gene_cluster_matrix_raw = standardize_cluster_matrix(gene_cluster_matrix_raw, gene_col=self.gene_col_name,
                                                             normalizations=self.normalizations)

        all_clusters = gene_cluster_matrix_raw.columns.tolist()
        all_genes = gene_cluster_matrix_raw[self.gene_col_name].tolist()

        gene_info_df = retrieve_gene_list_and_optional_annotations(
            gene_select_mode=self.gene_select_mode,
            fetch=fetch,
        )
        cluster_list,cluster_annotations= retrieve_cluster_list_and_optional_annotations(
            cluster_select_mode=self.cluster_select_mode,
            fetch=fetch,
            all_clusters=all_clusters,
        )

        # assert len(set(gene_info_df[MAGMA_GENE_COL]).intersection( set(all_genes)))>0
        assert set(cluster_list) <= set(all_clusters)


        gene_cluster_matrix_filtered = gene_cluster_matrix_raw.copy()
        # gene_cluster_matrix =gene_cluster_matrix.set_index(self.gene_col_name)
        gene_cluster_matrix_filtered = gene_cluster_matrix_filtered.merge(gene_info_df, left_on=self.gene_col_name, right_on=MAGMA_GENE_COL)#
        gene_cluster_matrix_filtered =gene_cluster_matrix_filtered.drop(columns=[self.gene_col_name, MAGMA_GENE_COL, MAGMA_P_COLUMN]).set_index(GENE_LABEL_COL)

        gene_cluster_matrix_filtered = gene_cluster_matrix_filtered.loc[:,gene_cluster_matrix_filtered.columns.isin(cluster_list)]

        vlim_params = get_vlim_params(cluster_matrix_filtered=gene_cluster_matrix_filtered,
                                      spec=self.vlimspec)

        logger.debug(f"Plotting cluster matrix with dimensions {gene_cluster_matrix_filtered.shape}")



        result = sns.clustermap(
            data=gene_cluster_matrix_filtered,
        # center = 0,
        #     cmap = "vlag",vmin=-2, vmax=2,
            **vlim_params,
            **self.plot_params
        )
        result.savefig(
            scratch_dir/"expression_clustermap.png",

        )
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(cls,
               spec_matrix_task: Task,
               asset_id: str,
               gene_select_mode: GeneSelectMode,
               cluster_select_mode: ClusterSelectMode,
               gene_col_name: str,
               normalizations: list[NormalizationStrategy],
               plot_params: dict,
               vlimspec: Vlimspec|None,
               ):
        if isinstance(gene_select_mode, TopGenesMode):
            src_meta = gene_select_mode.source_task.meta
            meta = GWASPlotDirectoryMeta(
                trait=src_meta.trait,
                project=src_meta.project,
                short_id=AssetId(asset_id)
            )
            return cls(
                meta=meta,
                specificity_matrix_task=spec_matrix_task,
                gene_select_mode=gene_select_mode,
                cluster_select_mode=cluster_select_mode,
                gene_col_name=gene_col_name,
                normalizations=normalizations,
                plot_params=plot_params,
                vlimspec=vlimspec
            )
        raise NotImplementedError("Not implemented for this gene select mode")

GENE_LABEL_COL = "gene_label"


def get_sorted_gene_df(
        fetch: Fetch,
        gene_select_mode: GeneSelectMode,
) -> pd.DataFrame:
    gene_source_asset = fetch(gene_select_mode.source_task.asset_id)
    assert isinstance(gene_source_asset, DirectoryAsset)
    gene_df = read_magma_gene_analysis_result(gene_source_asset.path)
    gene_df = gene_df.sort_values(by=MAGMA_P_COLUMN)
    return gene_df


def retrieve_gene_list_and_optional_annotations(
        gene_select_mode:GeneSelectMode,
        fetch: Fetch,
)-> pd.DataFrame:
    if isinstance(gene_select_mode, TopGenesMode):
        gene_df = get_sorted_gene_df(fetch=fetch, gene_select_mode=gene_select_mode)
        gene_df = gene_df.loc[gene_df[MAGMA_P_COLUMN] <= (gene_select_mode.p_value_thresh / len(gene_df))]
        sorted_filtered  =gene_df.sort_values(by=MAGMA_P_COLUMN).loc[:,[ MAGMA_GENE_COL, MAGMA_P_COLUMN, MAGMA_Z_COLUMN ]]
        gene_label_list = [f"{nm} (z={zscore})" for nm,p_val, zscore in zip(sorted_filtered[MAGMA_GENE_COL], sorted_filtered[MAGMA_P_COLUMN], sorted_filtered[MAGMA_Z_COLUMN])]
        sorted_filtered[GENE_LABEL_COL] = gene_label_list
        sorted_filtered=sorted_filtered.sort_values(by=MAGMA_P_COLUMN).iloc[:gene_select_mode.num_genes]
        return sorted_filtered.loc[:,[MAGMA_GENE_COL, GENE_LABEL_COL, MAGMA_P_COLUMN]]
    if isinstance(gene_select_mode, GeneNamesMode):
        pass



def retrieve_cluster_list_and_optional_annotations(
        cluster_select_mode:ClusterSelectMode,
        fetch: Fetch,
        all_clusters: list[str]
)-> tuple[list[str], list[str]|None]:
    if isinstance(cluster_select_mode, AllClustersMode):
        return all_clusters,None
    raise NotImplementedError("Not implemented")


def standardize_cluster_matrix(
        cluster_matrix: pd.DataFrame,
        gene_col: str,
        normalizations: list[NormalizationStrategy]
) -> pd.DataFrame:
    df =cluster_matrix.set_index(gene_col)
    for norm in normalizations:
        if isinstance(norm, GeneWiseSpecificityNormalization):
            for row in df.index:
                df.loc[row,:] = df.loc[row,:] / df.loc[row,:].sum()
        elif isinstance(norm, ClusterWiseZScoreNormalization):
            for col in df.columns:
                df[col] = (df[col] - df[col].mean())/(df[col].std())
    return df.reset_index()


def get_vlim_params(cluster_matrix_filtered: pd.DataFrame, spec: Vlimspec | None)-> dict:
    if spec is None:
        return {}
    if isinstance(spec, VlimForGeneWiseSpec):
        num_cols = cluster_matrix_filtered.shape[1]
        expected_frac= 1/num_cols
        return {
            "center": expected_frac,
            "vmin":0,
            "vmax":expected_frac*spec.max_multiple
        }



def _debugging(
        gene_select_mode: GeneSelectMode,
        fetch: Fetch,
        gene_cluster_matrix_raw : pd.DataFrame,
):
    gene_df = get_sorted_gene_df(fetch=fetch, gene_select_mode=gene_select_mode)
    print("yo")
