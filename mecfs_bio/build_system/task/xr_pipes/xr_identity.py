import xarray as xr

from mecfs_bio.build_system.task.xr_pipes.xr_data_pipe import (
    XRDataPipe,
)


class XRIdentityPipe(XRDataPipe):
    def process(self, ds: xr.Dataset) -> xr.Dataset:
        return ds
