from pathlib import Path

import numpy as np
import plotly.express as px
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.plot_meta import GWASPlotDirectoryMeta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.multiple_testing_table_task import (
    REJECT_NULL_LABEL,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.util.plotting.save_fig import write_plots_to_dir


@frozen
class SLDSCScatterPlotTask(Task):
    """
    Generate a plot from the resulting of applying S-LDSC to GWAS summary statistics
    using GTEx and Franke lab reference data.

    Intended to mimic the appearance of the plots in
    Finucane, Hilary K., et al. "Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types."
    Nature genetics 50.4 (2018): 621-629.
    """

    _meta: Meta
    df_source_task: Task

    @property
    def _source_id(self) -> AssetId:
        return self.df_source_task.asset_id

    @property
    def _source_meta(self) -> Meta:
        return self.df_source_task.meta

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.df_source_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        df_asset = fetch(self._source_id)
        df = (
            scan_dataframe_asset(df_asset, meta=self._source_meta).collect().to_pandas()
        )
        df = df.sort_values(by=["Category", "Name"])
        df["mlog10p"] = -np.log10(df["Coefficient_P_value"])
        df["tissue_or_cell"] = np.arange(len(df))
        df["Category"] = df["Category"].str.title()
        df["Category"] = df["Category"].str.replace("Cns", "CNS")
        df["marker_size"] = 4.5 * df[REJECT_NULL_LABEL] + 0.5
        fig = px.scatter(
            df,
            x="tissue_or_cell",
            y="mlog10p",
            color="Category",
            hover_name="Name",
            size="marker_size",
            hover_data={
                "Name": True,
                "Category": True,
                "mlog10p": True,
                "marker_size": False,
                "tissue_or_cell": False,
            },
            template="plotly_white",
            labels={
               "mlog10p":"$-log\u2081\u2080p$",
                "tissue_or_cell":"Tissue or Cell"
            }
        )
        fig = fig.update_xaxes(showticklabels=False)
        # fig= fig.update_layout(xaxis_title=None)
        figs = {}
        figs["sldsc_scatter"] = fig
        out_dir = scratch_dir / "scatter_plots"
        write_plots_to_dir(out_dir, figs, plotly_mathjax_mode="cdn")
        return DirectoryAsset(out_dir)

    @classmethod
    def create(cls, asset_id: str, source_task: Task):
        source_meta = source_task.meta
        assert isinstance(source_meta, ResultTableMeta)
        meta = GWASPlotDirectoryMeta(
            trait=source_meta.trait,
            project=source_meta.project,
            id=AssetId(asset_id),
        )
        return cls(
            meta=meta,
            df_source_task=source_task,
        )
