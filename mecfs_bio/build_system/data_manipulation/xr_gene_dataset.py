import pandas as pd
import xarray as xr

XR_GENE_DIMENSION = "GENE"
XR_TISSUE_DIMENSION = "TISSUE" # represents tissues/ cell-types, etc. used in MAGMA-like tools


XR_SPECIFICITY_MATRIX = "specificity_matrix"



def gene_spec_df_to_da(gene_tissue_df: pd.DataFrame, gene_col:str) -> xr.DataArray:
    return xr.DataArray(
        gene_tissue_df.set_index(gene_col),dims=(XR_GENE_DIMENSION, XR_TISSUE_DIMENSION),
    )

def gene_annotation_df_to_ds(gene_annotation_df: pd.DataFrame, gene_col: str|None)-> xr.Dataset:
    gene_annotation_df =gene_annotation_df.copy()
    if gene_col is not None:
        gene_annotation_df = gene_annotation_df.set_index(gene_col)
    gene_annotation_df.index.name = XR_GENE_DIMENSION
    return gene_annotation_df.to_xarray()

def tissue_annotation_df_to_ds(tissue_annotation_df:pd.DataFrame, tissue_col: str|None)-> xr.Dataset:
    tissue_annotation_df =tissue_annotation_df.copy()
    if tissue_col is not None:
        tissue_annotation_df  =tissue_annotation_df.set_index(tissue_col)
    tissue_annotation_df.index.name = XR_TISSUE_DIMENSION
    return tissue_annotation_df.to_xarray()



# def cluster_genes_by_tissue(
#     ds: xr.Dataset,
#     metric: _MetricKind = "correlation",
#     method : _LinkageMethod="average"
#     ) -> xr.Dataset:
#     X = ds[XR_SPECIFICITY_MATRIX].values
#     row_dist = pdist(X, metric=metric)
#     row_linkage = linkage(row_dist, method=method)
#     row_order = leaves_list(row_linkage)
#     ds = ds.isel({
#         XR_GENE_DIMENSION: row_order,
#     })
#     return ds
#
# def cluster_tissues_by_gene(
#         ds: xr.Dataset,
#         metric: _MetricKind = "correlation",
#         method : _LinkageMethod="average"
# ) -> xr.Dataset:
#     X = ds[XR_SPECIFICITY_MATRIX].values.T
#     col_dist = pdist(X, metric=metric)
#     col_linkage = linkage(col_dist, method=method)
#     col_order = leaves_list(col_linkage)
#     ds = ds.isel({
#         XR_TISSUE_DIMENSION: col_order,
#     })
#     return ds
#
