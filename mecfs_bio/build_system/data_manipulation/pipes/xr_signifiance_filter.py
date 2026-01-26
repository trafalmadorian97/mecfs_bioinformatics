import xarray as xr
from attrs import frozen

from mecfs_bio.build_system.data_manipulation.pipes.xr_data_pipe import XRDataPipe
from mecfs_bio.build_system.data_manipulation.xr_gene_dataset import XR_GENE_DIMENSION


@frozen
class XRSignifianceFilter(XRDataPipe):
    p_threshold: float
    p_da: str
    z_da:str
    count_dim: str = XR_GENE_DIMENSION

    def process(self, ds: xr.Dataset) -> xr.Dataset:
        count = ds.sizes[self.count_dim]
        adjusted_thresh = self.p_threshold/count
        ds= ds.where( (ds[self.p_da] <=adjusted_thresh) &  (ds[self.z_da]>=0), drop=True )
        return ds