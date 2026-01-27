import pandas as pd
import xarray as xr

from mecfs_bio.constants.xr_constants import XR_GENE_DIMENSION, XR_TISSUE_DIMENSION


def gene_spec_df_to_da(gene_tissue_df: pd.DataFrame, gene_col: str) -> xr.DataArray:
    return xr.DataArray(
        gene_tissue_df.set_index(gene_col),
        dims=(XR_GENE_DIMENSION, XR_TISSUE_DIMENSION),
    )


def gene_annotation_df_to_ds(
    gene_annotation_df: pd.DataFrame, gene_col: str | None
) -> xr.Dataset:
    gene_annotation_df = gene_annotation_df.copy()
    if gene_col is not None:
        gene_annotation_df = gene_annotation_df.set_index(gene_col)
    gene_annotation_df.index.name = XR_GENE_DIMENSION
    return gene_annotation_df.to_xarray()


def tissue_annotation_df_to_ds(
    tissue_annotation_df: pd.DataFrame, tissue_col: str | None
) -> xr.Dataset:
    tissue_annotation_df = tissue_annotation_df.copy()
    if tissue_col is not None:
        tissue_annotation_df = tissue_annotation_df.set_index(tissue_col)
    tissue_annotation_df.index.name = XR_TISSUE_DIMENSION
    return tissue_annotation_df.to_xarray()
