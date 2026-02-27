"""
Task to  plot the results of applying MAGMA using the human brain atlas as a reference.
Attempts to match the style of the graphs from
    Duncan, Laramie E., et al. "Mapping the cellular etiology of schizophrenia and complex brain phenotypes." Nature Neuroscience 28.2 (2025): 248-258.
"""

import math
from pathlib import Path, PurePath
from typing import Literal

import numpy as np
import pandas as pd
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
PlotMode = Literal["plotly_white", "plotly_dark"]

KEY_HBA_ANNOTATION_COLUMNS = [
    "Supercluster",
    "Class auto-annotation",
    "Neurotransmitter auto-annotation",
    "Neuropeptide auto-annotation",
    "Subtype auto-annotation",
    "Transferred MTG Label",
    "Top three regions",
    "Top Enriched Genes",
]


@frozen
class PlotSettings:
    plot_mode: PlotMode = "plotly_dark"


@frozen
class PlotMagmaBrainAtlasResultTask(Task):
    """
    Task to create a plot of the results of applying MAGMA using the human brain atlas data as a reference

    Attempts to match the style of the graphs from:
    Duncan, Laramie E., et al. "Mapping the cellular etiology of schizophrenia and complex brain phenotypes." Nature Neuroscience 28.2 (2025): 248-258.

    Note that "cluster annotation task" brings in some metadata, but this metadata is currently not used in the
    graph. In the future, we may use this metadata
    """

    result_table_task: Task
    cluster_annotation_task: Task
    _meta: Meta
    plot_mode: PlotMode = "plotly_dark"

    @property
    def annotation_id(self) -> AssetId:
        return self.cluster_annotation_task.asset_id

    @property
    def source_id(self) -> AssetId:
        return self.result_table_task.asset_id

    @property
    def annotation_meta(self) -> Meta:
        return self.cluster_annotation_task.meta

    @property
    def source_meta(self) -> Meta:
        return self.result_table_task.meta

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.result_table_task, self.cluster_annotation_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self.source_id)
        annotation_asset = fetch(self.annotation_id)
        table = (
            scan_dataframe_asset(asset=source_asset, meta=self.source_meta)
            .collect()
            .to_pandas()
        )
        annotation_table = (
            scan_dataframe_asset(asset=annotation_asset, meta=self.annotation_meta)
            .collect()
            .to_pandas()
        )
        table["MLOG10P"] = -np.log10(table["P"])
        table["CLUSTER"] = table["VARIABLE"].str.extract(r"Cluster([0-9]+)").astype(int)
        table = table.sort_values(by="CLUSTER")
        sorted_table = table.sort_values(by="P")
        top_cluster = sorted_table.iloc[0]
        top_cluster_label = get_condensed_hba_cluster_label(top_cluster)
        colormap = {
            key.replace(".", " "): val
            for key, val in DUNCAN_ET_AL_2025_COLORMAP.items()
        }
        plot = px.scatter(
            table,
            x="CLUSTER",
            y="MLOG10P",
            labels={
                "CLUSTER": "Cluster",
                "MLOG10P": "-log\u2081\u2080p",
            },
            color="Supercluster",
            color_discrete_map=colormap,
            template=self.plot_mode,
            hover_data={
                "Supercluster": True,
                "Class auto-annotation": True,
                "Neurotransmitter auto-annotation": True,
                "Neuropeptide auto-annotation": True,
                "Subtype auto-annotation": True,
                "Transferred MTG Label": True,
                "Top three regions": True,
                "Top Enriched Genes": True,
            },
            # size_max=20
        )
        plot = plot.update_traces(marker=dict(size=15))
        plot.update_layout(
            xaxis_title_font=dict(size=20),
            yaxis_title_font=dict(size=23),
        )

        sig_level = -math.log10(0.01 / len(table))
        line_color = "white" if self.plot_mode == "plotly_dark" else "black"
        plot = plot.add_hline(
            y=sig_level,
            line_color=line_color,
            line_dash="dot",
            opacity=0.5,
            annotation_text="Significance threshold (Bonferroni)",
            annotation_xshift=-10,
            annotation_font=dict(size=17),
        )
        plot = plot.add_annotation(
            x=float(top_cluster["CLUSTER"]),
            y=float(top_cluster["MLOG10P"]),
            text=top_cluster_label,
            font=dict(  # Customize font properties
                size=20,
                color=colormap[top_cluster["Supercluster"]],
            ),
            arrowhead=2,  # Style of the arrowhead
            arrowsize=1,  # Size of the arrowhead
            arrowwidth=2,  # Width of the arrow line
            arrowcolor=colormap[top_cluster["Supercluster"]],
            ay=-60,  # Y-component of the arrow tail offset (pixels from head)
            standoff=10,
        )
        plots = {BRAIN_ATLAS_PLOT_NAME: plot}
        write_plots_to_dir(scratch_dir, plots)
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        result_table_task: Task,
        cluster_annotation_task: Task,
        asset_id: str,
        plot_settings: PlotSettings = PlotSettings(),
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
            cluster_annotation_task=cluster_annotation_task,
            plot_mode=plot_settings.plot_mode,
        )


def add_cluster_column_to_metadata_df(metadata_df: pd.DataFrame) -> pd.DataFrame:
    # filter to get only clusters
    cluster_df = metadata_df.loc[
        metadata_df["cluster_annotation_term_set_name"] == "cluster"
    ]
    cluster_df["CLUSTER"] = (
        cluster_df["description"].str.extract("cluster ([0-9]+)").astype(int)
    )
    return cluster_df


"""
Colormap from Duncan et al. paper.
"""
DUNCAN_ET_AL_2025_COLORMAP = {
    # Non-neuronal / glial
    "Miscellaneous": "#8D4517",
    "Microglia": "#333333",
    "Vascular": "#4D4D4D",
    "Fibroblast": "#656565",
    "Oligodendrocyte.precursor": "#727272",
    "Committed.oligodendrocyte.precursor": "#7F7F7F",
    "Oligodendrocyte": "#8C8C8C",
    "Bergmann.glia": "#9A9A9A",
    "Astrocyte": "#A8A8A8",
    "Ependymal": "#BFBFBF",
    "Choroid.plexus": "#D9D9D9",
    # Excitatory neurons
    "Deep-layer.near-projecting": "#E41A1C",
    "Deep-layer.corticothalamic.and.6b": "#FF7F00",
    "Hippocampal.CA1-3": "#FDBF6F",
    "Upper-layer.intratelencephalic": "#FFFF33",
    "Deep-layer.intratelencephalic": "#B2DF8A",
    "Amygdala.excitatory": "#66C2A5",
    "Hippocampal.CA4": "#33A02C",
    "Hippocampal.dentate.gyrus": "#1B9E77",
    # Striatal / medium spiny
    "Medium.spiny.neuron": "#00CC99",
    "Eccentric.medium.spiny.neuron": "#1CE6B3",
    # Other neuronal types
    "Splatter": "#00E5FF",
    "MGE.interneuron": "#1F78B4",
    "LAMP5-LHX6.and.Chandelier": "#6A3D9A",
    "CGE.interneuron": "#0000FF",
    "Upper.rhombic.lip": "#5E3C99",
    "Cerebellar.inhibitory": "#7B3294",
    "Lower.rhombic.lip": "#8E0152",
    "Mammillary.body": "#E7298A",
    "Thalamic.excitatory": "#FF1493",
    "Midbrain-derived.inhibitory": "#E31A1C",
}


def get_condensed_hba_cluster_label(cluster_info: pd.Series) -> str:
    """
    Get a short label for a cluster than can be used to annotate points on the plot
    """
    result = str(cluster_info["Supercluster"])
    extra_info = []
    if cluster_info["Class auto-annotation"] not in ["0", "NEUR"]:
        extra_info.append(cluster_info["Class auto-annotation"])
    if cluster_info["Subtype auto-annotation"] != "0":
        extra_info.append(cluster_info["Subtype auto-annotation"])
    if cluster_info["Neurotransmitter auto-annotation"] != "0":
        extra_info.append(cluster_info["Neurotransmitter auto-annotation"])
    if (
        cluster_info["Neuropeptide auto-annotation"] != "0"
        and len(cluster_info["Neuropeptide auto-annotation"]) < 6
    ):
        extra_info.append(cluster_info["Neuropeptide auto-annotation"])
    if len(extra_info) > 0:
        result += f" ({', '.join(extra_info)})"
    return result
