from abc import ABC, abstractmethod

import xarray as xr


class XRDataPipe(ABC):
    """
    A transformation of an xarray datast
    """

    @abstractmethod
    def process(self, ds: xr.Dataset) -> xr.Dataset:
        pass
