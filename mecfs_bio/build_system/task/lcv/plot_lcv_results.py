from pathlib import Path

from attrs import frozen

from mecfs_bio.asset_generator.lcv_asset_generator import UPSTREAM_TRAIT_COL, DOWNSTREAM_TRAIT_COL
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.lcv.lcv_core import LCV_RHO_EST_COL, LCV_PVAL_ZERO_COL, LCV_MEAN_GCP_COL
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.xr_pipes.xr_data_pipe import XRDataPipe
from mecfs_bio.build_system.wf.base_wf import WF
import xarray as xr


XR_UPSTREAM_TRAIT_DIM = "upstream_trait"
XR_DOWNSTREAM_TRAIT_DIM = "downstream_trait"

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
    gcp_col: str= LCV_MEAN_GCP_COL
    df_pipe: DataProcessingPipe = IdentityPipe()
    @property
    def cols(self)-> list[str]:
        return [self.upstream_trait_col,
                self.downstream_trait_col,
                self.rho_col,
                self.p_col,
                self.gcp_col
                ]



def load_xr_corr_dataset(
        src: LCVSource,
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
            *src.cols,
        )
        .collect()
    )
    # num_pairs = _count_unique_pairs(df_nw_stacked, src.trait_1_col, src.trait_2_col)
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
    gcp_da = xr.DataArray(pivoted_gcp, dims=(XR_UPSTREAM_TRAIT_DIM   , XR_DOWNSTREAM_TRAIT_DIM))
    p_da = xr.DataArray(pivoted_p, dims=(XR_UPSTREAM_TRAIT_DIM   , XR_DOWNSTREAM_TRAIT_DIM))
    rho_da = xr.DataArray(pivoted_rho, dims=(XR_UPSTREAM_TRAIT_DIM   , XR_DOWNSTREAM_TRAIT_DIM))
    
    ds = xr.Dataset({XR_GENETIC_CORR_ARRAY: rg_da, XR_GENETIC_CORR_P_VALUE_ARRAY: p_da})
    # ds[NUM_PAIRS] = num_pairs
    return ds


@frozen
class LCVClustermapTask(Task):
    """
    Task to generate a heatmap of LCV REsults
    """

    meta: Meta
    xr_pipe: XRDataPipe
    source: LCVSource
    # genetic_corr_source: GeneticCorrSource
    # plot_options: GeneticCorrPlotMode
    save_mode: bool | str = "cdn"



    @property
    def deps(self) -> list["Task"]:
        return [self.source.task]


    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        # ds = load_xr_corr_dataset(
        #     src=self.genetic_corr_source,
        #     fetch=fetch,
        # )
        # ds = self.xr_pipe.process(ds)
        # fig = rg_plot(
        #     ds=ds,
        #     plot_mode=self.plot_options,
        # )
        # out_path = scratch_dirtch_dir / "result.html"
        fig.write_html(out_path, include_plotlyjs=self.save_mode)
        return FileAsset(out_path)
