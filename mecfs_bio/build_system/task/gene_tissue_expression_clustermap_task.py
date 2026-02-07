"""
Task to make a heatmap plot with genes as rows and tissue/cell types as columns
"""

from typing import Sequence

import seaborn as sns
import structlog

from mecfs_bio.build_system.data_manipulation.xr_data.pipes.xr_cluster import XRCluster
from mecfs_bio.build_system.data_manipulation.xr_data.pipes.xr_composite import (
    XRCompositePipe,
)
from mecfs_bio.build_system.data_manipulation.xr_data.pipes.xr_data_pipe import (
    XRDataPipe,
)
from mecfs_bio.build_system.data_manipulation.xr_data.pipes.xr_most_significant import (
    XRMostSignificant,
)
from mecfs_bio.build_system.data_manipulation.xr_data.pipes.xr_significance_filter import (
    XRSignificanceFilter,
)
from mecfs_bio.build_system.data_manipulation.xr_data.xr_gene_dataset_load import (
    GeneInfoSource,
    SpecificityMatrixSource,
    TissueInfoSource,
    load_xr_gene_tissue_dataset,
)
from mecfs_bio.build_system.data_manipulation.xr_data.xr_gene_specificity_heatmap import (
    GeneNormalizedPlotMode,
    XRHeatmapPlotSpec,
    xr_make_gene_spec_heatmap,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.plot_meta import GWASPlotDirectoryMeta
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.constants.magma_constants import MAGMA_P_COLUMN, MAGMA_Z_COLUMN
from mecfs_bio.constants.xr_constants import (
    XR_GENE_DIMENSION,
    XR_SPECIFICITY_MATRIX,
    XR_TISSUE_DIMENSION,
)
from mecfs_bio.util.plotting.save_fig import write_plots_to_dir

logger = structlog.get_logger()

from pathlib import Path

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class ExpressionMatrixClusterMapTaskV2(Task):
    """
    Task to make a heatmap plot with genes as rows and tissue/cell types as columns
    Goal is to show the specificity of genes for tissue/cell types

    """

    _meta: Meta
    specificity_matrix_source: SpecificityMatrixSource
    gene_info_sources: Sequence[GeneInfoSource]
    tissue_info_sources: Sequence[TissueInfoSource]
    xr_pipe: XRDataPipe
    plot_spec: XRHeatmapPlotSpec

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        result = [self.specificity_matrix_source.task]
        for gi in self.gene_info_sources:
            result.append(gi.task)
        for ti in self.tissue_info_sources:
            result.append(ti.task)
        return result

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        ds = load_xr_gene_tissue_dataset(
            gt_specificity_matrix_source=self.specificity_matrix_source,
            gene_info_sources=self.gene_info_sources,
            tissue_info_sources=self.tissue_info_sources,
            fetch=fetch,
        )
        ds = self.xr_pipe.process(ds)
        fig = xr_make_gene_spec_heatmap(
            ds=ds,
            plot_spec=self.plot_spec,
        )
        write_plots_to_dir(scratch_dir, {"gene_tissue_heatmap": fig})
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create_standard_gene_magma_heatmap(
        cls,
        asset_id: str,
        extracted_magma_gene_results_source: GeneInfoSource,
        gene_specificity_matrix_source: SpecificityMatrixSource,
        gene_thesaurus_source: GeneInfoSource,
        num_genes_to_keep: int = 300,
    ):
        """
        Make a heatmap in which genes are selected according to their MAGMA significance.
        Useful for delving deeper into MAGMA results .
        """
        source_meta = extracted_magma_gene_results_source.task.meta
        assert isinstance(source_meta, ResultTableMeta)
        meta = GWASPlotDirectoryMeta(
            id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
        )
        seaborn_palette = sns.color_palette("vlag", n_colors=100)
        plotly_colorscale = [
            f"rgb({int(c[0] * 255)}, {int(c[1] * 255)}, {int(c[2] * 255)})"
            for c in seaborn_palette
        ]
        return cls(
            meta=meta,
            specificity_matrix_source=gene_specificity_matrix_source,
            gene_info_sources=[
                extracted_magma_gene_results_source,
                gene_thesaurus_source,
            ],
            tissue_info_sources=[],
            xr_pipe=XRCompositePipe(
                [
                    XRCluster(
                        array_name=XR_SPECIFICITY_MATRIX,
                        dim=XR_TISSUE_DIMENSION,
                    ),
                    XRSignificanceFilter(
                        p_threshold=0.01,
                        p_da=MAGMA_P_COLUMN,
                        z_da=MAGMA_Z_COLUMN,
                    ),
                    XRMostSignificant(
                        ordering_da=MAGMA_P_COLUMN,
                        num_to_keep=num_genes_to_keep,
                    ),
                    XRCluster(
                        array_name=XR_SPECIFICITY_MATRIX,
                        dim=XR_GENE_DIMENSION,
                    ),
                ]
            ),
            plot_spec=XRHeatmapPlotSpec(
                GeneNormalizedPlotMode(max_multiple=5),
                y_array="Gene name",
                color_scale=plotly_colorscale,
            ),
        )
