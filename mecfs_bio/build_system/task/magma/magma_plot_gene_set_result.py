from pathlib import Path, PurePath

import numpy as np
import pandas as pd
import plotly.express as px
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.plot_meta import GWASPlotDirectoryMeta
from mecfs_bio.build_system.meta.processed_gwas_data_directory_meta import (
    ProcessedGwasDataDirectoryMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.magma.magma_gene_set_analysis_task import (
    GENE_SET_ANALYSIS_OUTPUT_STEM_NAME,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import GWASLAB_MLOG10P_COL
from mecfs_bio.util.plotting.save_fig import write_plots_to_dir

logger = structlog.get_logger()

MAGMA_GENE_SET_PLOT_NAME = "magma_gene_set_plot"


@frozen
class MAGMAPlotGeneSetResult(Task):
    _meta: Meta
    gene_set_analysis_task: Task
    number_of_bars: int = 20

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def _source_id(self) -> AssetId:
        return self.gene_set_analysis_task.asset_id

    @property
    def deps(self) -> list["Task"]:
        return [self.gene_set_analysis_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        gene_set_dir_asset = fetch(self._source_id)
        assert isinstance(gene_set_dir_asset, DirectoryAsset)
        result_path = gene_set_dir_asset.path / str(
            GENE_SET_ANALYSIS_OUTPUT_STEM_NAME + ".gsa.out"
        )
        df = pd.read_csv(result_path, comment="#", sep=r"\s+")
        num_tests = len(df)
        df[GWASLAB_MLOG10P_COL] = -np.log10(df["P"])
        logger.debug(f"Found {num_tests} tests")
        thresh = 0.05 / num_tests
        logger.debug(f"Bonferroni thresh: {thresh}")
        trunc_df = df.sort_values(by=["P"], ascending=True).iloc[: self.number_of_bars]
        trunc_df["Significance"] = float("nan")
        trunc_df["Significance"] = "Significant"
        trunc_df["Significance"] = trunc_df["Significance"].where(
            trunc_df["P"] <= thresh, "Not Significant"
        )
        plot = px.bar(
            trunc_df,
            x="FULL_NAME",
            y=GWASLAB_MLOG10P_COL,
            labels={
                "FULL_NAME": "Variable",
                GWASLAB_MLOG10P_COL: "-log\u2081\u2080p",
            },
            color="Significance",
            color_discrete_map={
                "Significant": px.colors.qualitative.Plotly[1],
                "Not Significant": px.colors.qualitative.Plotly[0],
            },
        )
        plot.add_hline(
            -np.log10(thresh),
            line_dash="dot",
            annotation_text="Significance threshold (Bonferroni)",
            annotation_xshift=-10,
        )
        plots = {MAGMA_GENE_SET_PLOT_NAME: plot}
        write_plots_to_dir(scratch_dir, plots)
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls, gene_set_analysis_task: Task, asset_id: str, number_of_bars: int = 20
    ):
        source_meta = gene_set_analysis_task.meta
        assert isinstance(source_meta, ProcessedGwasDataDirectoryMeta)
        meta = GWASPlotDirectoryMeta(
            trait=source_meta.trait,
            project=source_meta.project,
            short_id=AssetId(asset_id),
            sub_dir=PurePath("analysis/magma_plots"),
        )
        return cls(
            meta,
            gene_set_analysis_task=gene_set_analysis_task,
            number_of_bars=number_of_bars,
        )
