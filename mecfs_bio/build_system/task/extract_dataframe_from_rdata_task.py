from pathlib import Path, PurePath
from typing import Sequence

import pandas as pd
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec, DataFrameParquetFormat
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import ReferenceFileMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

import pyreadr
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import (
    importr,
)
import rpy2.robjects as ro

@frozen
class ExtractDataFrameFromRDataTask(Task):
    """
    RData files can bundle together many R objects.  This is a Task to extract a single dataframe from such a file.
    """
    _meta: Meta
    rdata_file_task: Task
    r_dataframe_name:str
    r_package_list:Sequence[str]

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.rdata_file_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        for package in self.r_package_list:
            importr(package)
        conv = ro.default_converter + pandas2ri.converter + numpy2ri.converter
        rdata_asset = fetch(self.rdata_file_task.asset_id)
        assert isinstance(rdata_asset, FileAsset)
        pth = rdata_asset.path
        robjects.r['load'](str(pth))
        r_dataframe = robjects.r[self.r_dataframe_name]
        with localconverter(conv):
            py_dataframe:pd.DataFrame = ro.conversion.get_conversion().rpy2py(r_dataframe)
        out_path= scratch_dir/"df.parquet"
        py_dataframe.to_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(cls,
               asset_id:str,
        rdata_file_task: Task,
        r_dataframe_name: str,
        r_package_list: Sequence[str]
    ):
        source_meta = rdata_file_task.meta
        if isinstance(source_meta, ReferenceFileMeta):
            meta= ReferenceFileMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                extension=".parquet",
                read_spec=DataFrameReadSpec(
                    DataFrameParquetFormat()
                )
            )
        else:
            raise ValueError(f"Unknown meta: {source_meta}")

        return cls(
            meta=meta,
            rdata_file_task=rdata_file_task,
            r_dataframe_name=r_dataframe_name,
            r_package_list=r_package_list
        )


