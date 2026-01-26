from typing import Literal

import xarray as xr
from attrs import frozen

from scipy.cluster.hierarchy import linkage,  leaves_list
from scipy.spatial.distance import pdist
from mecfs_bio.build_system.data_manipulation.pipes.xr_data_pipe import XRDataPipe
import structlog

Linkage = Literal["average"]
Metric = Literal["correlation"]
logger = structlog.get_logger()

@frozen
class XRCluster(XRDataPipe):
    array_name: str
    dim:str
    metric: Metric= "correlation"
    method : Linkage="average"

    def process(self, ds: xr.Dataset) -> xr.Dataset:
        array = ds[self.array_name]
        array = array.transpose(self.dim,...)
        X=array.values

        logger.debug(f"Clustering X with X.shape={X.shape}")
        row_dist = pdist(X, metric=self.metric)
        row_linkage = linkage(row_dist, method=self.method)
        row_order = leaves_list(row_linkage)
        ds = ds.isel({
            self.dim: row_order,
        })
        return ds

