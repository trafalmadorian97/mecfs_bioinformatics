from pathlib import Path
from zipfile import ZipFile

from attrs import frozen
from tqdm import tqdm

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class ExtractAllZipsInDir(Task):
    _meta: Meta
    source_directory_task: Task

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.source_directory_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self.source_directory_task.asset_id)
        assert isinstance(source_asset, DirectoryAsset)
        target_dir = scratch_dir / "extracted"
        target_dir.mkdir(parents=True, exist_ok=True)
        for pth in tqdm(sorted(source_asset.path.glob("*.zip"))):
            with ZipFile(pth, "r") as zip_file:
                zip_file.extractall(
                    target_dir,
                )
        return DirectoryAsset(target_dir)

    @classmethod
    def create(cls, source_directory_task: Task, asset_id: str):
        source_meta = source_directory_task.meta
        if isinstance(source_meta, ReferenceDataDirectoryMeta):
            meta = ReferenceDataDirectoryMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=source_meta.sub_folder,
                id=AssetId(asset_id),
            )
        else:
            raise ValueError(f"Unknown meta {source_meta} ")
        return cls(
            source_directory_task=source_directory_task,
            meta=meta,
        )
