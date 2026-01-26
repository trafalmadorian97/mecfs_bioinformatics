import xarray as xr
from attrs import frozen

from mecfs_bio.build_system.data_manipulation.pipes.xr_data_pipe import XRDataPipe


@frozen
class XRMostSignificant(XRDataPipe):
    ordering_da: str
    num_to_keep:int

    def process(self, ds: xr.Dataset) -> xr.Dataset:
        target_dim = ds[self.ordering_da].dims[0]
        ds = ds.sortby(self.ordering_da).isel({target_dim: slice(None, self.num_to_keep )})
        return ds
