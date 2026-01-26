import xarray as xr
from abc import ABC, abstractmethod


class XRDataPipe(ABC):
    @abstractmethod
    def process(self, ds: xr.Dataset) -> xr.Dataset:
        pass