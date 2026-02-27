import xarray as xr

from mecfs_bio.build_system.data_manipulation.xr_data.pipes.xr_data_pipe import (
    XRDataPipe,
)


class XRIdentityPipe(XRDataPipe):
    def process(self, ds: xr.Dataset) -> xr.Dataset:
        return ds
