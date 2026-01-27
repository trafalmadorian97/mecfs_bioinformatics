from typing import Sequence

import structlog
import xarray as xr
from attrs import frozen

from mecfs_bio.build_system.data_manipulation.xr_data.pipes.xr_data_pipe import (
    XRDataPipe,
)

logger = structlog.get_logger()


@frozen
class XRCompositePipe(XRDataPipe):
    """
    Combine multiple XRDataPipes
    """

    pipes: Sequence[XRDataPipe]

    def process(self, ds: xr.Dataset) -> xr.Dataset:
        for pipe in self.pipes:
            logger.debug(f"Running pipe {pipe}")
            ds = pipe.process(ds)
        return ds
