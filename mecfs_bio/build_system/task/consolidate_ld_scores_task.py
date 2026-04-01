"""
Task to read LD scores in the standard format defined by the authors of LD score regression,
and write them out as a parquet file
"""

from pathlib import Path

import narwhals
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class ConsolidateLDScoresTask(Task):
    """
    Task to read LD scores in the standard format defined by the authors of LD score regression,
    and write them out as a parquet file
    """

    _meta: Meta
    extracted_ld_scores_task: Task

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.extracted_ld_scores_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        asset = fetch(self.extracted_ld_scores_task.asset_id)
        assert isinstance(asset, DirectoryAsset)
        frames = []
        for ld_file in sorted(asset.path.glob("*.gz")):
            frames.append(
                scan_dataframe(
                    ld_file,
                    DataFrameReadSpec(DataFrameTextFormat(separator="\t")),
                )
            )
        result = narwhals.concat(frames, how="vertical").sort(by=["CHR", "BP"])
        out_path = scratch_dir / "out.parquet"
        result.sink_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(cls, asset_id: str, extracted_ld_score_task: Task):
        source_meta = extracted_ld_score_task.meta
        assert isinstance(source_meta, ReferenceDataDirectoryMeta)
        meta = ReferenceFileMeta(
            group=source_meta.group,
            sub_group=source_meta.sub_group,
            sub_folder=source_meta.sub_folder,
            id=AssetId(asset_id),
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            extension=".parquet",
        )
        return cls(meta=meta, extracted_ld_scores_task=extracted_ld_score_task)
