import math
from pathlib import Path, PurePath

import numpy as np
import plotly.express as px
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.plot_meta import GWASPlotDirectoryMeta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import (
    scan_dataframe_asset,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.util.plotting.save_fig import write_plots_to_dir

BRAIN_ATLAS_PLOT_NAME = "human_brain_atlas_plot"


@frozen
class PlotMagmaBrainAtlasResultTask(Task):
    result_table_task: Task
    _meta: Meta

    @property
    def source_id(self) -> AssetId:
        return self.result_table_task.asset_id

    @property
    def source_meta(self) -> Meta:
        return self.result_table_task.meta

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.result_table_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self.source_id)
        table = (
            scan_dataframe_asset(asset=source_asset, meta=self.source_meta)
            .collect()
            .to_pandas()
        )
        table["MLOG10P"] = -np.log10(table["P"])

        table["CLUSTER"] = table["VARIABLE"].str.extract(r"Cluster([0-9]+)").astype(int)
        table = table.sort_values(by="CLUSTER")
        plot = px.scatter(
            table,
            x="CLUSTER",
            y="MLOG10P",
            color="Supercluster",
            hover_data={
                "Supercluster": True,
                "Class auto-annotation": True,
                "Neurotransmitter auto-annotation": True,
                "Neuropeptide auto-annotation": True,
                "Subtype auto-annotation": True,
                "Transferred MTG Label": True,
                "Top three regions": True,
            },
        )
        sig_level = -math.log10(0.01 / len(table))
        plot = plot.add_hline(
            y=sig_level,
            line_color="black",
            # annotation_text=f"Significance Level: {sig_level}",
            line_dash="dash",
        )
        plots = {BRAIN_ATLAS_PLOT_NAME: plot}
        write_plots_to_dir(scratch_dir, plots)
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        result_table_task: Task,
        asset_id: str,
    ):
        source_meta = result_table_task.meta
        assert isinstance(source_meta, ResultTableMeta)
        meta = GWASPlotDirectoryMeta(
            trait=source_meta.trait,
            project=source_meta.project,
            short_id=AssetId(asset_id),
            sub_dir=PurePath("analysis/magma_plots"),
        )
        return cls(meta=meta, result_table_task=result_table_task)
