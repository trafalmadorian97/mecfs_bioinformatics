import math
from pathlib import Path, PurePath

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import textalloc as ta
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
from mecfs_bio.build_system.task.magma.magma_forward_stepwise_select_task import (
    RETAINED_CLUSTERS_COLUMN,
)
from mecfs_bio.build_system.task.magma.plot_magma_brain_atlas_result import (
    DUNCAN_ET_AL_2025_COLORMAP,
    get_condensed_hba_cluster_label,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.util.plotting.save_fig import write_plots_to_dir


@frozen
class HBAIndepPlotOptions:
    annotation_text_size: int = 7


@frozen
class MAGMAPlotBrainAtlasResultWithStepwiseLabels(Task):
    """
    Create a plot of the results of applying MAGMA to a GWAS using HBA reference data

    - Since the gene expression patterns of brain cell types are highly correlated, it can be difficult
     to distinguish the independent significant signals on a typical plot of MAGMA results.  For this reason, this plotting a function
     annotates a subset of the cell types that are relatively independent significant: i.e. the p value of one cell type
     does not decline too much when conditioning on the others.
    """

    result_table_task: Task
    stepwise_cluster_list_task: Task
    _meta: Meta
    plot_options: HBAIndepPlotOptions = HBAIndepPlotOptions()

    @property
    def stepwise_clusters_id(
        self,
    ) -> AssetId:
        return self.stepwise_cluster_list_task.asset_id

    @property
    def stepwise_met(self) -> Meta:
        return self.stepwise_cluster_list_task.meta

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
        return [self.result_table_task, self.stepwise_cluster_list_task]

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

        sig_level = -math.log10(0.01 / len(table))

        stepwise_clusters_asset = fetch(self.stepwise_clusters_id)
        independent_clusters = (
            scan_dataframe_asset(asset=stepwise_clusters_asset, meta=self.stepwise_met)
            .collect()
            .to_pandas()[RETAINED_CLUSTERS_COLUMN]
            .tolist()
        )

        colormap = {
            key.replace(".", " "): val
            for key, val in DUNCAN_ET_AL_2025_COLORMAP.items()
        }
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16, 10))
        sns.scatterplot(
            table,
            x="CLUSTER",
            y="MLOG10P",
            hue="Supercluster",
            palette=colormap,
            ax=ax,
            s=145,
        )
        plt.subplots_adjust(
            left=0.05,
            right=0.8,
        )
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1.0, 1), fontsize="small")
        ax.axhline(y=sig_level, color="black", linestyle="--", linewidth=4, zorder=-100)
        ax.set_xlabel("CLUSTER", fontsize="x-large")
        ax.set_ylabel("MLOG10P", fontsize="x-large")

        if len(independent_clusters) > 0:
            x_coords = []
            y_cords = []
            texts = []
            for item in independent_clusters:
                x_coords.append(table.loc[table["VARIABLE"] == item]["CLUSTER"].item())
                y_cords.append(table.loc[table["VARIABLE"] == item]["MLOG10P"].item())
                texts.append(
                    get_condensed_hba_cluster_label(
                        table.loc[table["VARIABLE"] == item].iloc[0, :],
                    )
                )
            ta.allocate(
                ax,
                x_coords,
                y_cords,
                texts,
                x_scatter=table["CLUSTER"].tolist(),
                y_scatter=table["MLOG10P"].tolist(),
                x_lines=[[table["CLUSTER"].max(), table["CLUSTER"].min()]],
                y_lines=[[sig_level, sig_level]],
                linewidth=1,
                avoid_crossing_label_lines=True,
                avoid_label_lines_overlap=True,
                linecolor="black",
                textsize=self.plot_options.annotation_text_size,
            )

        figs = {"hba_magma_fig": fig}
        write_plots_to_dir(scratch_dir, figs)

        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        result_table_task: Task,
        asset_id: str,
        stepwise_cluster_list_task: Task,
        plot_options: HBAIndepPlotOptions = HBAIndepPlotOptions(),
    ):
        source_meta = result_table_task.meta
        assert isinstance(source_meta, ResultTableMeta)
        meta = GWASPlotDirectoryMeta(
            trait=source_meta.trait,
            project=source_meta.project,
            id=AssetId(asset_id),
            sub_dir=PurePath("analysis/magma_plots"),
        )
        return cls(
            meta=meta,
            result_table_task=result_table_task,
            stepwise_cluster_list_task=stepwise_cluster_list_task,
            plot_options=plot_options,
        )
