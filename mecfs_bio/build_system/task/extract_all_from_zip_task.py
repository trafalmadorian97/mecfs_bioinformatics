from pathlib import Path, PurePath
from zipfile import ZipFile

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
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
class ExtractAllFromZipTask(Task):
    _meta: Meta
    source_file_task: Task

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def _source_meta(self) -> Meta:
        return self.source_file_task.meta

    @property
    def _source_id(self) -> AssetId:
        return self.source_file_task.asset_id

    @property
    def deps(self) -> list["Task"]:
        return [self.source_file_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self._source_id)
        assert isinstance(source_asset, FileAsset)
        target_dir = scratch_dir / "extracted"
        target_dir.mkdir(parents=True, exist_ok=True)
        with ZipFile(source_asset.path, "r") as zip_file:
            zip_file.extractall(
                target_dir,
            )
        return DirectoryAsset(target_dir)

    @classmethod
    def create_from_zipped_reference_file(
        cls, source_task: Task, asset_id: str
    ) -> Task:
        src_meta = source_task.meta
        assert isinstance(src_meta, ReferenceFileMeta)
        meta = ReferenceDataDirectoryMeta(
            group=src_meta.group,
            sub_group=src_meta.sub_group,
            sub_folder=PurePath("extracted"),
            id=AssetId(asset_id),
            dirname=src_meta.filename,
        )
        return cls(
            meta=meta,
            source_file_task=source_task,
        )
