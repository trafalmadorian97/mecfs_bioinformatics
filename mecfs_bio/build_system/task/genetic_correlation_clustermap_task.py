"""
Task to create a heatmap plot illustrating the genetic correlation structure of a collection of traits

Based on fig 2 from

Bulik-Sullivan, Brendan, et al. "An atlas of genetic correlations across human diseases and traits." Nature genetics 47.11 (2015): 1236-1241.

"""

from pathlib import Path

import narwhals
import numpy as np
import plotly.graph_objs as go
import structlog
import xarray as xr
from attrs import frozen
from plotly.graph_objs import Figure

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.data_manipulation.xr_data.pipes.xr_cluster import XRCluster
from mecfs_bio.build_system.data_manipulation.xr_data.pipes.xr_composite import (
    XRCompositePipe,
)
from mecfs_bio.build_system.data_manipulation.xr_data.pipes.xr_data_pipe import (
    XRDataPipe,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.plot_file_meta import GWASPlotFileMeta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()
XR_GENETIC_CORR_ARRAY = "genetic_correlation"
XR_GENETIC_CORR_P_VALUE_ARRAY = "genetic_correlation_p_value"
XR_TRAIT_1_DIM = "trait_1"
XR_TRAIT_2_DIM = "trait_2"
NUM_PAIRS = "num_pairs"


@frozen
class GeneticCorrSource:
    """
    Describe a dataframe from which to load genetic correlation data
    """

    task: Task
    trait_1_col: str = "p1"
    trait_2_col: str = "p2"
    rg_col: str = "rg"
    p_col: str = "p"
    df_pipe: DataProcessingPipe = IdentityPipe()


def _count_unique_pairs(df: narwhals.DataFrame, col1: str, col2: str) -> int:
    pairs = set(
        [
            frozenset([df[i, col1], df[i, col2]])
            for i in range(len(df))
            if df[i, col1] != df[i, col2]
        ]
    )
    return len(pairs)


def fill_out_corr_df(
    df_nw: narwhals.DataFrame,
    src: GeneticCorrSource,
) -> narwhals.DataFrame:
    """
    Full out the upper triangle and diagonal of the correlation matrix
    """
    df_2_nw = df_nw.with_columns(
        narwhals.col(src.trait_1_col).alias(src.trait_2_col),
        narwhals.col(src.trait_2_col).alias(
            src.trait_1_col
        ),  # fill out upper triangle of correlation matrix
    )
    df_nw = narwhals.concat([df_nw, df_2_nw], how="vertical").unique(keep="first")
    df_diag = (
        df_nw.select(src.trait_1_col)
        .unique()
        .with_columns(
            narwhals.col(src.trait_1_col).alias(src.trait_2_col),
            narwhals.lit(1.0).alias(src.rg_col),
            narwhals.lit(float("nan")).alias(src.p_col),
        )
    )  # fill out diagonal of correlation matrix
    df_nw = narwhals.concat([df_nw, df_diag], how="vertical").unique()
    return df_nw


def load_xr_corr_dataset(
    src: GeneticCorrSource,
    fetch: Fetch,
) -> xr.Dataset:
    """
    Retrieve the correlation data.  Return in the form of an xarray dataset.
    """
    asset = fetch(src.task.asset_id)
    df_nw = (
        src.df_pipe.process(
            scan_dataframe_asset(
                asset,
                meta=src.task.meta,
            )
        )
        .select(
            src.trait_1_col,
            src.trait_2_col,
            src.rg_col,
            src.p_col,
        )
        .collect()
    )
    df_nw_stacked = fill_out_corr_df(df_nw, src)
    num_pairs = _count_unique_pairs(df_nw_stacked, src.trait_1_col, src.trait_2_col)
    df = df_nw_stacked.to_pandas()

    pivoted_rg = df.pivot(
        index=src.trait_1_col,
        columns=src.trait_2_col,
        values=src.rg_col,
    )
    pivoted_p = df.pivot(
        index=src.trait_1_col,
        columns=src.trait_2_col,
        values=src.p_col,
    )
    rg_da = xr.DataArray(pivoted_rg, dims=(XR_TRAIT_1_DIM, XR_TRAIT_2_DIM))
    p_da = xr.DataArray(pivoted_p, dims=(XR_TRAIT_1_DIM, XR_TRAIT_2_DIM))
    ds = xr.Dataset({XR_GENETIC_CORR_ARRAY: rg_da, XR_GENETIC_CORR_P_VALUE_ARRAY: p_da})
    ds[NUM_PAIRS] = num_pairs
    return ds


@frozen
class BonferoniSig:
    alpha: float = 0.05


SigMode = BonferoniSig


@frozen
class RGWithAsterix:
    sig_mode: SigMode = BonferoniSig()
    color_scale: str = "RdBu_r"


@frozen
class RGHideNonSig:
    pass


GeneticCorrPlotMode = RGWithAsterix | RGHideNonSig


def get_sig(
    p_value_matrix: np.ndarray,
    num_pairs: int,
    sig_mode: SigMode,
) -> np.ndarray:
    """
    Get a binary array indicating which elements of the correlation matrix are significance
    """
    if isinstance(sig_mode, BonferoniSig):
        thresh = sig_mode.alpha / num_pairs
        logger.debug(
            f"Bonferoni significance threshold when alpha={sig_mode.alpha} and there are {num_pairs} tests is {thresh}"
        )
        return p_value_matrix <= thresh
    raise ValueError(f"Invalid mode {sig_mode}")


def rg_plot(ds: xr.Dataset, plot_mode: GeneticCorrPlotMode) -> Figure:
    """
    Produce a plotly heatmap figure showing genetic correlation
    """
    if isinstance(plot_mode, RGWithAsterix):
        corr_df = ds[XR_GENETIC_CORR_ARRAY].to_pandas()
        sig = get_sig(
            p_value_matrix=ds[XR_GENETIC_CORR_P_VALUE_ARRAY].values,
            num_pairs=ds[NUM_PAIRS].values.item(),
            sig_mode=plot_mode.sig_mode,
        )
        asterix_matrix = np.where(sig, "★", "")

        fig = go.Figure(
            data=go.Heatmap(
                z=corr_df,
                x=corr_df.columns,
                y=corr_df.index,
                text=asterix_matrix,
                texttemplate="%{text}",
                textfont={"size": 20, "color": "black"},
                colorscale=plot_mode.color_scale,
                zmin=-1,
                zmax=1,
                customdata=ds[XR_GENETIC_CORR_P_VALUE_ARRAY].values.tolist(),
                hovertemplate="Trait 1: %{y}<br>"
                + "Trait 2: %{x}<br>"
                + "Genetic Correlation: %{z}<br>"
                + "p value: %{customdata}<br>"
                + "<extra></extra>",  # Removes the default "trace 0" info
                showscale=True,
                hovertemplatefallback="None",
            )
        )
        fig.update_layout(
            xaxis=dict(side="top"),
            title="Genetic Correlations",
            xaxis_title="Trait 2",
            yaxis_title="Trait 1",
        )
        fig.update_yaxes(autorange="reversed")  # want the origin in top left corner
        return fig
    raise NotImplementedError()


@frozen
class GeneticCorrelationClustermapTask(Task):
    """
    Task to generate a heatmap of genetic correlation
    """

    _meta: Meta
    xr_pipe: XRDataPipe
    genetic_corr_source: GeneticCorrSource
    plot_options: GeneticCorrPlotMode
    save_mode: bool | str = "cdn"

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.genetic_corr_source.task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        ds = load_xr_corr_dataset(
            src=self.genetic_corr_source,
            fetch=fetch,
        )
        ds = self.xr_pipe.process(ds)
        fig = rg_plot(
            ds=ds,
            plot_mode=self.plot_options,
        )
        out_path = scratch_dir / "result.html"
        fig.write_html(out_path, include_plotlyjs=self.save_mode)
        return FileAsset(out_path)

    @classmethod
    def create_std_with_clustering(
        cls,
        asset_id: str,
        genetic_corr_source: GeneticCorrSource,
        plot_options: GeneticCorrPlotMode,
    ):
        src_meta = genetic_corr_source.task.meta
        assert isinstance(src_meta, ResultTableMeta)
        meta = GWASPlotFileMeta(
            trait=src_meta.trait,
            project=src_meta.project,
            extension=".html",
            id=AssetId(asset_id),
        )
        xr_pipe = XRCompositePipe(
            [
                XRCluster(
                    array_name=XR_GENETIC_CORR_ARRAY,
                    dim=XR_TRAIT_1_DIM,
                    metric="euclidean",
                ),
                XRCluster(
                    array_name=XR_GENETIC_CORR_ARRAY,
                    dim=XR_TRAIT_2_DIM,
                    metric="euclidean",
                ),
            ]
        )
        return cls(
            meta=meta,
            xr_pipe=xr_pipe,
            genetic_corr_source=genetic_corr_source,
            plot_options=plot_options,
        )
