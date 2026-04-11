"""
Task to create a heatmap plot illustrating the results of Latent Causal Variable
(LCV) analysis for a collection of upstream / downstream trait pairs.

The primary value plotted is the posterior mean of the Genetic Causality
Proportion (GCP). Significantly non-zero entries (under a Bonferroni correction
to the LCV GCP=0 p-values) are marked with an asterisk.

Partly implemented by asking Claude to mimic the logic of  `genetic_correlation_clustermap_task.py`.


"""

from pathlib import Path

import numpy as np
import plotly.graph_objs as go
import structlog
import xarray as xr
from attrs import frozen
from plotly.graph_objs import Figure

from mecfs_bio.asset_generator.lcv_asset_generator import (
    DOWNSTREAM_TRAIT_COL,
    UPSTREAM_TRAIT_COL,
)
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.plot_file_meta import GWASPlotFileMeta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.lcv.lcv_core import (
    LCV_MEAN_GCP_COL,
    LCV_PVAL_ZERO_COL,
    LCV_RHO_EST_COL,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.xr_pipes.xr_cluster import XRCluster
from mecfs_bio.build_system.task.xr_pipes.xr_composite import XRCompositePipe
from mecfs_bio.build_system.task.xr_pipes.xr_data_pipe import XRDataPipe
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()

XR_UPSTREAM_TRAIT_DIM = "upstream_trait"
XR_DOWNSTREAM_TRAIT_DIM = "downstream_trait"
XR_GCP_ARRAY = "gcp"
XR_LCV_P_VALUE_ARRAY = "lcv_p_value"
XR_LCV_RHO_ARRAY = "lcv_rho"
NUM_PAIRS = "num_pairs"


@frozen
class LCVSource:
    """
    Describe a dataframe from which to load LCV data
    """

    task: Task
    upstream_trait_col: str = UPSTREAM_TRAIT_COL
    downstream_trait_col: str = DOWNSTREAM_TRAIT_COL
    rho_col: str = LCV_RHO_EST_COL
    p_col: str = LCV_PVAL_ZERO_COL
    gcp_col: str = LCV_MEAN_GCP_COL
    df_pipe: DataProcessingPipe = IdentityPipe()

    @property
    def cols(self) -> list[str]:
        return [
            self.upstream_trait_col,
            self.downstream_trait_col,
            self.rho_col,
            self.p_col,
            self.gcp_col,
        ]


def load_xr_lcv_dataset(
    src: LCVSource,
    fetch: Fetch,
) -> xr.Dataset:
    """
    Retrieve the LCV results.  Return in the form of an xarray dataset.
    """
    asset = fetch(src.task.asset_id)
    df_nw = (
        src.df_pipe.process(
            scan_dataframe_asset(
                asset,
                meta=src.task.meta,
            )
        )
        .select(*src.cols)
        .collect()
    )
    num_pairs = len(df_nw)
    df = df_nw.to_pandas()

    pivoted_gcp = df.pivot(
        index=src.upstream_trait_col,
        columns=src.downstream_trait_col,
        values=src.gcp_col,
    )
    pivoted_rho = df.pivot(
        index=src.upstream_trait_col,
        columns=src.downstream_trait_col,
        values=src.rho_col,
    )
    pivoted_p = df.pivot(
        index=src.upstream_trait_col,
        columns=src.downstream_trait_col,
        values=src.p_col,
    )
    gcp_da = xr.DataArray(
        pivoted_gcp, dims=(XR_UPSTREAM_TRAIT_DIM, XR_DOWNSTREAM_TRAIT_DIM)
    )
    p_da = xr.DataArray(
        pivoted_p, dims=(XR_UPSTREAM_TRAIT_DIM, XR_DOWNSTREAM_TRAIT_DIM)
    )
    rho_da = xr.DataArray(
        pivoted_rho, dims=(XR_UPSTREAM_TRAIT_DIM, XR_DOWNSTREAM_TRAIT_DIM)
    )

    ds = xr.Dataset(
        {
            XR_GCP_ARRAY: gcp_da,
            XR_LCV_P_VALUE_ARRAY: p_da,
            XR_LCV_RHO_ARRAY: rho_da,
        }
    )
    ds[NUM_PAIRS] = num_pairs
    return ds


@frozen
class BonferoniSig:
    alpha: float = 0.05


SigMode = BonferoniSig


@frozen
class GCPWithAsterisk:
    sig_mode: SigMode = BonferoniSig()
    color_scale: str = "RdBu_r"


LCVPlotMode = GCPWithAsterisk


def get_sig(
    p_value_matrix: np.ndarray,
    num_pairs: int,
    sig_mode: SigMode,
) -> np.ndarray:
    """
    Get a binary array indicating which elements of the LCV result matrix are
    significant under the given significance scheme.
    """
    if isinstance(sig_mode, BonferoniSig):
        thresh = sig_mode.alpha / num_pairs
        logger.debug(
            f"Bonferoni significance threshold when alpha={sig_mode.alpha} "
            f"and there are {num_pairs} tests is {thresh}"
        )
        return p_value_matrix <= thresh
    raise ValueError(f"Invalid mode {sig_mode}")


def gcp_plot(ds: xr.Dataset, plot_mode: LCVPlotMode) -> Figure:
    """
    Produce a plotly heatmap figure showing LCV GCP estimates for each
    upstream / downstream trait pair.
    """
    if isinstance(plot_mode, GCPWithAsterisk):
        gcp_df = ds[XR_GCP_ARRAY].to_pandas()
        p_values = ds[XR_LCV_P_VALUE_ARRAY].values
        rho_values = ds[XR_LCV_RHO_ARRAY].values
        sig = get_sig(
            p_value_matrix=p_values,
            num_pairs=ds[NUM_PAIRS].values.item(),
            sig_mode=plot_mode.sig_mode,
        )
        asterisk_matrix = np.where(sig, "★", "")
        customdata = np.stack([p_values, rho_values], axis=-1)

        fig = go.Figure(
            data=go.Heatmap(
                z=gcp_df,
                x=gcp_df.columns,
                y=gcp_df.index,
                text=asterisk_matrix,
                texttemplate="%{text}",
                textfont={"size": 20, "color": "black"},
                colorscale=plot_mode.color_scale,
                zmin=-1,
                zmax=1,
                customdata=customdata,
                hovertemplate="Upstream Trait: %{y}<br>"
                + "Downstream Trait: %{x}<br>"
                + "GCP: %{z}<br>"
                + "Genetic correlation (rho): %{customdata[1]}<br>"
                + "p value (GCP=0): %{customdata[0]}<br>"
                + "<extra></extra>",
                showscale=True,
                hovertemplatefallback="None",
            )
        )
        fig.update_layout(
            xaxis=dict(side="top"),
        )
        fig.update_yaxes(autorange="reversed")  # want the origin in top left corner
        return fig
    raise NotImplementedError()


@frozen
class LCVClustermapTask(Task):
    """
    Task to generate a heatmap of LCV results
    """

    meta: Meta
    xr_pipe: XRDataPipe
    source: LCVSource
    plot_options: LCVPlotMode
    save_mode: bool | str = "cdn"

    @property
    def deps(self) -> list["Task"]:
        return [self.source.task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        ds = load_xr_lcv_dataset(
            src=self.source,
            fetch=fetch,
        )
        ds = self.xr_pipe.process(ds)
        fig = gcp_plot(
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
        source: LCVSource,
        plot_options: LCVPlotMode,
    ):
        src_meta = source.task.meta
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
                    array_name=XR_GCP_ARRAY,
                    dim=XR_UPSTREAM_TRAIT_DIM,
                    metric="euclidean",
                ),
                XRCluster(
                    array_name=XR_GCP_ARRAY,
                    dim=XR_DOWNSTREAM_TRAIT_DIM,
                    metric="euclidean",
                ),
            ]
        )
        return cls(
            meta=meta,
            xr_pipe=xr_pipe,
            source=source,
            plot_options=plot_options,
        )
