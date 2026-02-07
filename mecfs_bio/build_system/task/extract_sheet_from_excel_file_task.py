"""
Task to extract sheet from excel file.
"""

from pathlib import Path, PurePath
from typing import Mapping

import pandas as pd
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    CSVOutFormat,
    OutFormat,
    ParquetOutFormat,
    get_extension_and_read_spec_from_format,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class ExtractSheetFromExelFileTask(Task):
    """
    Task for extracting a sheet from an excel file
    Useful for working with supplementary material to research papers, which is often in excel format
    """

    _meta: Meta
    excel_file_task: Task
    sheet_name: str
    out_format: OutFormat
    skiprows: int | None = None
    col_type_mapping: Mapping[str, type] | None = None
    post_pipe: DataProcessingPipe = IdentityPipe()

    @property
    def source_asset_id(self) -> AssetId:
        return self.excel_file_task.asset_id

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.excel_file_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        excel_asset = fetch(self.source_asset_id)
        assert isinstance(excel_asset, FileAsset)
        out_path = scratch_dir / self.asset_id
        df = pd.read_excel(
            excel_asset.path, sheet_name=self.sheet_name, skiprows=self.skiprows
        )
        if self.col_type_mapping is not None:
            for col in self.col_type_mapping.keys():
                df[col] = df[col].astype(self.col_type_mapping[col])

        df = self.post_pipe.process_pandas(df)
        if isinstance(self.out_format, CSVOutFormat):
            df.to_csv(out_path, index=False, sep=self.out_format.sep)
        elif isinstance(self.out_format, ParquetOutFormat):
            df.to_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        source_excel_file_task: Task,
        sheet_name: str,
        out_format: OutFormat,
        skiprows: int | None = None,
        col_type_mapping: Mapping[str, type] | None = None,
        post_pipe: DataProcessingPipe = IdentityPipe(),
    ):
        extension, read_spec = get_extension_and_read_spec_from_format(
            out_format=out_format
        )
        source_meta = source_excel_file_task.meta

        if isinstance(source_meta, ReferenceFileMeta):
            meta = ReferenceFileMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                filename=asset_id,
                extension=extension,
                read_spec=read_spec,
            )
        else:
            raise ValueError("Unknown source meta")

        return cls(
            meta=meta,
            excel_file_task=source_excel_file_task,
            sheet_name=sheet_name,
            out_format=out_format,
            skiprows=skiprows,
            col_type_mapping=col_type_mapping,
            post_pipe=post_pipe,
        )
