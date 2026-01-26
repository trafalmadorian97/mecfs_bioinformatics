import xarray as xr
from attrs import frozen


@frozen
class MostSignificantSelection:
    p_threshold: float
    p_col: str
    z_col:str


def select(
        ds: xr.Dataset,
)-> xr.Dataset:
    pass
