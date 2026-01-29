from typing import Literal

import structlog

from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.processed_gwas_data_directory_meta import ProcessedGwasDataDirectoryMeta

logger = structlog.get_logger()
import tarfile
import tempfile
from pathlib import Path, PurePath

from attrs import frozen

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

ReadMode = Literal["r", "r:gz"]


@frozen
class ExtractTarGzipTask(Task):
    """
    Task to extract the contents of a (possibly gzipped) tar file to a target directory
    Set subdir_name to extract only the contents of one subfolder within the tar file


    read_mode: use r to only untar, and not ungzip.
    """

    _meta: Meta
    _source_file_task: Task
    _subdir_name: str | None
    _read_mode: ReadMode = "r:gz"

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self._source_file_task]

    @property
    def _source_asset_id(self) -> AssetId:
        return self._source_file_task.asset_id

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> DirectoryAsset:
        source_asset = fetch(self._source_asset_id)
        assert isinstance(source_asset, FileAsset)
        src_path = source_asset.path

        logger.debug(f"Extracting from tar/gzip file : {self._source_asset_id}...")
        with tarfile.open(src_path, self._read_mode) as tar_object:
            if self._subdir_name is None:
                tar_object.extractall(scratch_dir)
            else:
                with tempfile.TemporaryDirectory() as tmpdir_name:
                    tmpdir_path = Path(tmpdir_name)
                    for member in tar_object.getmembers():
                        if member.name.startswith(self._subdir_name):
                            tar_object.extract(member=member, path=tmpdir_path)
                    (tmpdir_path / self._subdir_name).rename(scratch_dir)
        logger.debug("Extraction complete.")
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        source_task: Task,
        sub_folder: PurePath = PurePath("extracted"),
        sub_folder_name_inside_tar: str | None = None,
        read_mode: ReadMode = "r:gz",
    ) -> "ExtractTarGzipTask":
        source_meta = source_task.meta
        if isinstance(source_meta, ReferenceFileMeta):
            return cls(
                meta=ReferenceDataDirectoryMeta(
                    group=source_meta.group,
                    sub_group=source_meta.sub_group,
                    sub_folder=sub_folder,
                    asset_id=AssetId(asset_id),
                ),
                source_file_task=source_task,
                subdir_name=sub_folder_name_inside_tar,
                read_mode=read_mode,
            )
        if isinstance(source_meta, GWASSummaryDataFileMeta):
            return cls(
                meta= ProcessedGwasDataDirectoryMeta(
                    short_id=AssetId(asset_id),
                    trait=source_meta.trait,
                    project=source_meta.project,
                    sub_dir=sub_folder,
                ),
                source_file_task=source_task,
                subdir_name=sub_folder_name_inside_tar,
                read_mode=read_mode,
            )
        raise NotImplementedError(f"Handler for meta {source_meta} not implemented")

