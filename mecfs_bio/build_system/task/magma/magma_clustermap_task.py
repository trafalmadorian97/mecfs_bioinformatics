import plotly.express as px
import sys
from typing import List, Sequence
import structlog
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.drop_col_pipe import DropColPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.util.plotting.save_fig import write_plots_to_dir

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


# todo:
# Try rewriting to work with xarray datasets
#  advantage: can carry zscores, p values, alternate gene names, specificity matrix all in same dataset
#  Slicing genes will slice everything associated with genes
#  cleaner than having to repeatedly merge.

@frozen
class TopGenesMode:
    num_genes: int
    source_task: MagmaGeneAnalysisTask
    p_value_thresh:float

@frozen
class GeneNamesMode:
    gene_names: List[str]


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
class ExtraSeabornPlotOptions:
    hide_dendogram: bool
    hide_gene_labels:bool
    rotate_cluster_labels:bool



# class PlotlyDivergingFromMidpointColor:
#     color_scale: str



# PlotlyColorOptions = PlotlyDivergingFromMidpointColor

@frozen
class ExtraPlotlyImshowPlotOptions:
    pass
    # color_options : PlotlyColorOptions

PlotOptions = ExtraSeabornPlotOptions| ExtraPlotlyImshowPlotOptions


@frozen
class AuxGeneNamingDFSpec:
    aux_info_task:Task

@frozen
class GeneLabelingSpec:
    gene_name_merge_col:str|None
    name_to_use_col: str|None



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
    specificity_matrix_pipe: DataProcessingPipe
    gene_select_mode:GeneSelectMode
    cluster_select_mode: ClusterSelectMode
    gene_col_name: str
    normalizations: list[NormalizationStrategy]
    extra_plot_options: PlotOptions
    plot_params: dict
    gene_labeling_spec: GeneLabelingSpec
    vlimspec: Vlimspec|None
    aux_gene_naming_df_spec: AuxGeneNamingDFSpec | None = None



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
        if self.aux_gene_naming_df_spec is not None:
            result.append(self.aux_gene_naming_df_spec.aux_info_task)
        return result

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        gene_cluster_raw_asset = fetch(self.specificity_matrix_task.asset_id)

        gene_cluster_matrix_raw = self.specificity_matrix_pipe.process(scan_dataframe_asset(gene_cluster_raw_asset,meta=self.specificity_matrix_task.meta,)).collect().to_pandas()
        # _debugging(
        #     gene_select_mode=self.gene_select_mode,
        #     fetch=fetch,
        #     gene_cluster_matrix_raw=gene_cluster_matrix_raw,
        # )

        gene_cluster_matrix_raw = standardize_cluster_matrix(gene_cluster_matrix_raw, gene_col=self.gene_col_name,
                                                             normalizations=self.normalizations)
        gene_cluster_matrix_raw = drop_constant_genes(gene_cluster_matrix_raw, gene_col=self.gene_col_name)

        gene_cluster_matrix_raw = cluster_gene_dataframe(gene_cluster_matrix_raw, gene_col=self.gene_col_name,
                                                         cluster_genes=False,
                                                         cluster_clusters=True)

        all_clusters = gene_cluster_matrix_raw.columns.tolist()
        all_genes = gene_cluster_matrix_raw[self.gene_col_name].tolist()

        gene_info_df = retrieve_gene_list_and_optional_annotations(
            gene_select_mode=self.gene_select_mode,
            fetch=fetch,
            gene_subset=all_genes,
        )
        cluster_list,cluster_annotations= retrieve_cluster_list_and_optional_annotations(
            cluster_select_mode=self.cluster_select_mode,
            fetch=fetch,
            all_clusters=all_clusters,
        )

        # assert len(set(gene_info_df[MAGMA_GENE_COL]).intersection( set(all_genes)))>0
        assert set(cluster_list) <= set(all_clusters)


        gene_cluster_matrix_filtered = gene_cluster_matrix_raw.copy()
        gene_cluster_matrix_filtered = gene_cluster_matrix_filtered.merge(gene_info_df, left_on=self.gene_col_name, right_on=MAGMA_GENE_COL)#


        aux_gene_naming_df = load_aux_gene_naming_df(spec=self.aux_gene_naming_df_spec, fetch=fetch)

        gene_cluster_matrix_filtered = prepare_cluster_matrix_with_label_col(
            gene_cluster_matrix_filtered=gene_cluster_matrix_filtered, gene_col_name=self.gene_col_name,
            aux_naming_df=aux_gene_naming_df,
            gene_labeling=self.gene_labeling_spec
        )

        gene_cluster_matrix_filtered = cluster_gene_dataframe(
            gene_cluster_matrix=gene_cluster_matrix_filtered,
            gene_col=None,
            cluster_genes=True,
            cluster_clusters=False
        )


        gene_cluster_matrix_filtered = gene_cluster_matrix_filtered.loc[:,gene_cluster_matrix_filtered.columns.isin(cluster_list)]
        logger.debug(f"Plotting cluster matrix with dimensions {gene_cluster_matrix_filtered.shape}")
        write_fig(
            gene_cluster_matrix_filtered=gene_cluster_matrix_filtered,
            scratch_dir=scratch_dir,
            plot_params=self.plot_params,
            extra_plot_options=self.extra_plot_options,
            vlimspec=self.vlimspec,
        )





        return DirectoryAsset(scratch_dir)



    @classmethod
    def create(cls,
               spec_matrix_task: Task,
               spec_matrix_pipe: DataProcessingPipe,
               asset_id: str,
               gene_select_mode: GeneSelectMode,
               cluster_select_mode: ClusterSelectMode,
               gene_col_name: str,
               normalizations: list[NormalizationStrategy],
               plot_params: dict,
               extra_plot_options: PlotOptions,
               vlimspec: Vlimspec|None,
               gene_labeling_spec: GeneLabelingSpec,
               aug_gene_naming_df_spec: AuxGeneNamingDFSpec|None
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
                extra_plot_options=extra_plot_options,
                vlimspec=vlimspec,
                gene_labeling_spec=gene_labeling_spec,
                aux_gene_naming_df_spec=aug_gene_naming_df_spec,
                specificity_matrix_pipe=spec_matrix_pipe

            )
        raise NotImplementedError("Not implemented for this gene select mode")

    @classmethod
    def create_standard_seaborn(cls,
                                spec_matrix_task: Task,
                                top_genes_task: MagmaGeneAnalysisTask,
                                gene_thesaurus_task:Task,
                                cols_to_drop: Sequence[str]= ("Average",),
                                ):
        return cls.create(
    spec_matrix_task=spec_matrix_task,
    asset_id="pgc2022_sch_magma_clustermap_plotly",
    gene_select_mode=TopGenesMode(
        num_genes=300,
        source_task=top_genes_task,
        p_value_thresh=0.01
    ),
    cluster_select_mode=AllClustersMode(),
    gene_col_name="Gene",
    normalizations=[
        GeneWiseSpecificityNormalization(),
    ],
    plot_params={
             "cmap" : "vlag",
        "figsize":(20,11),
        "metric":'correlation',
            "dendrogram_ratio":(0.1,0.1)

    },
    vlimspec=VlimForGeneWiseSpec(max_multiple=5),
extra_plot_options=ExtraPlotlyImshowPlotOptions(
),
    gene_labeling_spec=GeneLabelingSpec("Gene stable ID","Gene name"),
    aug_gene_naming_df_spec=AuxGeneNamingDFSpec(gene_thesaurus_task),
    spec_matrix_pipe=DropColPipe(
        # []
        cols_to_drop=list(cols_to_drop),
    )

)

GENE_LABEL_COL = "gene_label"
GENE_AUX_INFO_COL = "gene_aux_info"


def write_fig(
              gene_cluster_matrix_filtered: pd.DataFrame,
              scratch_dir: Path,
             plot_params: dict,
            extra_plot_options: PlotOptions,
        vlimspec: Vlimspec|None,
              ) -> None:
    if isinstance(extra_plot_options, ExtraSeabornPlotOptions):
        vlim_params = get_vlim_seaborn_params(cluster_matrix_filtered=gene_cluster_matrix_filtered,
                                              spec=vlimspec)
        if extra_plot_options.hide_gene_labels:
            plot_params["yticklabels"]=False
        result = sns.clustermap(
            data=gene_cluster_matrix_filtered,
            **vlim_params,
            **plot_params
        )
        if extra_plot_options.hide_dendogram:
            result.ax_row_dendrogram.set_visible(False)
            result.ax_col_dendrogram.set_visible(False)
        if extra_plot_options.hide_gene_labels:
            result.ax_heatmap.set_ylabel("")
        if extra_plot_options.rotate_cluster_labels:
            plt.setp(
                result.ax_heatmap.get_xticklabels(),
                rotation=45,
                horizontalalignment='right'  # Aligns the end of the text with the tick
            )

        result.savefig(
            scratch_dir / "expression_clustermap.png",
        )
        return
    if isinstance(extra_plot_options, ExtraPlotlyImshowPlotOptions):
        vlim_params = get_vlim_plotly_params(
            cluster_matrix_filtered=gene_cluster_matrix_filtered,
            spec=vlimspec
        )
        fig = px.imshow(
            gene_cluster_matrix_filtered,
            color_continuous_scale="Tropic",
            # width=1800,
            # height=1100,
            aspect="auto",
            **vlim_params,
        )
        write_plots_to_dir(
            scratch_dir,
            {
                "expression_heatmap": fig
            }
        )
        return



def load_aux_gene_naming_df(
       spec: AuxGeneNamingDFSpec | None ,
        fetch: Fetch
)-> pd.DataFrame|None:
    if spec is None:
        return None
    asset = fetch(spec.aux_info_task.asset_id)
    df = scan_dataframe_asset(asset, meta=spec.aux_info_task.meta).collect().to_pandas()
    return df

def prepare_cluster_matrix_with_label_col(gene_cluster_matrix_filtered: pd.DataFrame,
                                          gene_col_name:str,
                                          aux_naming_df: pd.DataFrame | None,
                                          gene_labeling: GeneLabelingSpec
                                          ) -> pd.DataFrame:
    gene_cluster_matrix_filtered =gene_cluster_matrix_filtered.copy()
    if aux_naming_df  is None or gene_labeling.name_to_use_col is None:
        gene_cluster_matrix_filtered[GENE_LABEL_COL] = gene_cluster_matrix_filtered[MAGMA_GENE_COL] + \
                                                   gene_cluster_matrix_filtered[GENE_AUX_INFO_COL]

        gene_cluster_matrix_filtered = gene_cluster_matrix_filtered.drop(
            columns=[gene_col_name, MAGMA_GENE_COL, MAGMA_P_COLUMN, GENE_AUX_INFO_COL]).set_index(GENE_LABEL_COL)
        return gene_cluster_matrix_filtered

    aux_naming_df  = aux_naming_df.loc[:,[gene_labeling.name_to_use_col, gene_labeling.gene_name_merge_col]]
    aux_naming_df =aux_naming_df.drop_duplicates(subset=[gene_labeling.gene_name_merge_col])
    gene_cluster_matrix_filtered = gene_cluster_matrix_filtered.merge(aux_naming_df, how="left", left_on=gene_col_name, right_on=gene_labeling.gene_name_merge_col)

    gene_cluster_matrix_filtered[GENE_LABEL_COL] = gene_cluster_matrix_filtered[gene_labeling.name_to_use_col] + \
                                                   gene_cluster_matrix_filtered[GENE_AUX_INFO_COL]

    gene_cluster_matrix_filtered = gene_cluster_matrix_filtered.drop(
        columns=[gene_col_name, MAGMA_GENE_COL, MAGMA_P_COLUMN, GENE_AUX_INFO_COL, gene_labeling.name_to_use_col, gene_labeling.gene_name_merge_col]).set_index(GENE_LABEL_COL)

    return gene_cluster_matrix_filtered

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
        gene_subset: list[str]| None,
        fetch: Fetch,
)-> pd.DataFrame:
    if isinstance(gene_select_mode, TopGenesMode):
        gene_df = get_sorted_gene_df(fetch=fetch, gene_select_mode=gene_select_mode)
        if gene_subset is not None:
            gene_df = gene_df.loc[gene_df[MAGMA_GENE_COL].isin(gene_subset)]
        gene_df = gene_df.loc[(gene_df[MAGMA_P_COLUMN] <= (gene_select_mode.p_value_thresh / len(gene_df))) & gene_df[MAGMA_Z_COLUMN]>0]
        sorted_filtered  =gene_df.sort_values(by=MAGMA_P_COLUMN).loc[:,[ MAGMA_GENE_COL, MAGMA_P_COLUMN, MAGMA_Z_COLUMN ]]
        # gene_label_list = [f"{nm} (z={zscore})" for nm,p_val, zscore in zip(sorted_filtered[MAGMA_GENE_COL], sorted_filtered[MAGMA_P_COLUMN], sorted_filtered[MAGMA_Z_COLUMN])]
        gene_aux_list = [f"(z={zscore})" for nm,p_val, zscore in zip(sorted_filtered[MAGMA_GENE_COL], sorted_filtered[MAGMA_P_COLUMN], sorted_filtered[MAGMA_Z_COLUMN])]
        sorted_filtered[GENE_AUX_INFO_COL] = gene_aux_list
        sorted_filtered=sorted_filtered.sort_values(by=MAGMA_P_COLUMN).iloc[:gene_select_mode.num_genes]
        return sorted_filtered.loc[:,[MAGMA_GENE_COL, GENE_AUX_INFO_COL, MAGMA_P_COLUMN]]
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

def drop_constant_genes(
        cluster_matrix:pd.DataFrame,
        gene_col: str
) -> pd.DataFrame:

    cluster_matrix =cluster_matrix.set_index(gene_col)
    constant = cluster_matrix.std(axis=1) <=sys.float_info.epsilon
    cluster_matrix =cluster_matrix.loc[~constant]
    return cluster_matrix.reset_index()




def get_vlim_seaborn_params(cluster_matrix_filtered: pd.DataFrame, spec: Vlimspec | None)-> dict:
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

def get_vlim_plotly_params(cluster_matrix_filtered: pd.DataFrame, spec: Vlimspec | None)-> dict:
    if spec is None:
        return {}
    if isinstance(spec, VlimForGeneWiseSpec):
        num_cols = cluster_matrix_filtered.shape[1]
        expected_frac= 1/num_cols
        zmin = expected_frac-spec.max_multiple*expected_frac
        zmax = expected_frac + spec.max_multiple*expected_frac
        return {
            "zmin":zmin,
            "zmax":zmax,
        }


def _debugging(
        gene_select_mode: GeneSelectMode,
        fetch: Fetch,
        gene_cluster_matrix_raw : pd.DataFrame,
):
    gene_df = get_sorted_gene_df(fetch=fetch, gene_select_mode=gene_select_mode)
    print("yo")




def cluster_gene_dataframe(
        gene_cluster_matrix: pd.DataFrame,
        gene_col:str|None,
        cluster_genes: bool,
        cluster_clusters: bool,
)-> pd.DataFrame:
    if gene_col is not None:
        gene_cluster_matrix = gene_cluster_matrix.set_index(gene_col)
    clustered = cluster_dataframe(gene_cluster_matrix,
                                  cluster_rows=cluster_genes,
                                  cluster_cols=cluster_clusters)
    if gene_col is not None:
        clustered = clustered.reset_index()
    return clustered

def cluster_dataframe(
df: pd.DataFrame,
cluster_rows:bool,
cluster_cols: bool,
*,
row_metric: str = "correlation",
col_metric: str = "correlation",
method: str = "average",
) -> pd.DataFrame:
    """
    Hierarchically cluster rows and columns of a DataFrame and return
    a reordered copy.
    From Chatgpt
    """

    X = df.values


    if cluster_rows:
        logger.debug(f"clustering rows of dataframe of shape {df.shape}")
        row_dist = pdist(X, metric=row_metric)
        row_linkage = linkage(row_dist, method=method)
        row_order = leaves_list(row_linkage)
        df = df.iloc[row_order]


    if cluster_cols:

        logger.debug(f"clustering cols of dataframe of shape {df.shape}")
        col_dist = pdist(X.T, metric=col_metric)
        col_linkage = linkage(col_dist, method=method)
        col_order = leaves_list(col_linkage)
        df=df.iloc[:, col_order]



    logger.debug(f"Done clustering")


    return df