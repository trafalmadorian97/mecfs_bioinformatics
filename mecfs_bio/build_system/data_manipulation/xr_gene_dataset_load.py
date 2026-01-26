import xarray as xr
from typing import Sequence

from attrs import frozen

from mecfs_bio.build_system.data_manipulation.xr_gene_dataset import gene_spec_df_to_da, XR_SPECIFICITY_MATRIX, \
    gene_annotation_df_to_ds, tissue_annotation_df_to_ds
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class SpecificityMatrixSource:
    task: Task
    pipe: DataProcessingPipe
    gene_col:str

@frozen
class GeneInfoSource:
    task: Task
    pipe: DataProcessingPipe
    gene_col:str|None

@frozen
class TissueInfoSource:
    task: Task
    pipe: DataProcessingPipe
    tissue_col:str|None


def load_xr_gene_tissue_dataset(
        gt_specificity_matrix_source:SpecificityMatrixSource,
        gene_info_sources: Sequence[GeneInfoSource],
        tissue_info_sources: Sequence[TissueInfoSource],
        fetch: Fetch,
    )-> xr.Dataset:
    spec_asset =fetch(gt_specificity_matrix_source.task.asset_id)
    spec_df = gt_specificity_matrix_source.pipe.process( scan_dataframe_asset(spec_asset,
                                   meta=gt_specificity_matrix_source.task.meta,)).collect().to_pandas()

    spec_da=gene_spec_df_to_da(
        gene_tissue_df=spec_df,
        gene_col=gt_specificity_matrix_source.gene_col,
    )
    ds = xr.Dataset(
        {
        XR_SPECIFICITY_MATRIX: spec_da
        }
    )
    for gene_info_source in gene_info_sources:
        gene_info_asset = fetch(gene_info_source.task.asset_id)
        gene_info_df = gene_info_source.pipe.process(scan_dataframe_asset(gene_info_asset,meta=gene_info_source.task.meta,)).collect().to_pandas()
        gene_info_df = gene_info_df.drop_duplicates(subset=[gene_info_source.gene_col])
        gene_info_ds =gene_annotation_df_to_ds(
            gene_annotation_df=gene_info_df,
            gene_col=gene_info_source.gene_col,
        )
        ds = xr.merge([ds, gene_info_ds], join="left")

    for tissue_info_source in tissue_info_sources:
        tissue_info_asset = fetch(tissue_info_source.task.asset_id)
        tissue_info_df =tissue_info_source.pipe.process(
            scan_dataframe_asset(tissue_info_asset,meta=tissue_info_source.task.meta,)
        ).collect().to_pandas()
        tissue_info_ds = tissue_annotation_df_to_ds(
            tissue_annotation_df=tissue_info_df,
            tissue_col=tissue_info_source.tissue_col,
        )
        ds = xr.merge([ds, tissue_info_ds], join="left")

    return ds