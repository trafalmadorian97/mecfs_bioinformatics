"""
Create an effect measure plot using zepid
"""

from pathlib import Path

import matplotlib.pyplot as plt
from attrs import frozen
from zepid.graphics import EffectMeasurePlot

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.plot_meta import GWASPlotDirectoryMeta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class EffectMeasurePlotConfig:
    y_label_col: str
    y_label: str
    effect_size_col: str
    effect_size_label: str
    se_col: str
    ref_line_center: float
    t_adjust: float
    title: str | None = None
    figsize: tuple[float, float] = (7, 3)


@frozen
class PlotMREffectMeasure(Task):
    """
    Task to create an effect measure plot
    Useful for example to show the results of Mendelian Randomization

    See:
    https://zepid.readthedocs.io/en/latest/Reference/Graphics.html?highlight=effectmeasureplot
    """

    _meta: Meta
    source_df_task: Task
    config: EffectMeasurePlotConfig
    pre_pipe: DataProcessingPipe = IdentityPipe()

    @property
    def source_id(self) -> AssetId:
        return self.source_df_task.asset_id

    @property
    def source_meta(self) -> Meta:
        return self.source_df_task.meta

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.source_df_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self.source_id)
        df = (
            self.pre_pipe.process(
                scan_dataframe_asset(source_asset, meta=self.source_meta)
            )
            .collect()
            .to_pandas()
        )
        df = df.sort_values(by=self.config.effect_size_col)
        mean = df[self.config.effect_size_col]
        se = df[self.config.se_col]
        lower = mean - 1.95 * se
        upper = mean + 1.95 * se
        p = EffectMeasurePlot(
            label=df[self.config.y_label_col].tolist(),
            effect_measure=mean.tolist(),
            lcl=lower.tolist(),
            ucl=upper.tolist(),
        )
        p.labels(
            center=self.config.ref_line_center,
            effectmeasure=self.config.effect_size_label,
        )
        ax = p.plot(
            figsize=self.config.figsize,
            t_adjuster=self.config.t_adjust,
        )
        if self.config.title is not None:
            plt.title(self.config.title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(True)
        ax.spines["left"].set_visible(False)
        plt.savefig(scratch_dir / "effect_measure_plot.png", bbox_inches="tight")
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        source_df_task: Task,
        config: EffectMeasurePlotConfig,
        pre_pipe: DataProcessingPipe = IdentityPipe(),
    ):
        source_meta = source_df_task.meta
        assert isinstance(source_meta, ResultTableMeta)
        meta = GWASPlotDirectoryMeta(
            trait=source_meta.trait,
            project=source_meta.project,
            short_id=AssetId(asset_id),
        )
        return cls(
            source_df_task=source_df_task,
            config=config,
            meta=meta,
            pre_pipe=pre_pipe,
        )
