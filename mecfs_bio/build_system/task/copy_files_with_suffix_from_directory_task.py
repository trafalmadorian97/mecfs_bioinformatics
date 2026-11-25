"""
Copy every file with a given suffix (for example .svg) from the top level of a DirectoryAsset
into a new plot directory.

Useful when a task writes figures whose filenames are only known at build time, such as one
figure per locus, alongside other outputs: the copy holds just the figures in one format, so
they can be exported to the docs figure directory as a unit. A source with no matching files
gives an empty directory.
"""

import shutil
from pathlib import Path
from typing import Self

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.plot_meta import GWASPlotDirectoryMeta
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF


@frozen(slots=True)
class CopyFilesWithSuffixFromDirectoryTask(Task):
    """
    source_directory_task: a Task producing a DirectoryAsset.
    suffix: the file suffix to copy, including the leading dot (for example .svg).
    """

    meta: GWASPlotDirectoryMeta
    source_directory_task: Task
    suffix: str

    def __attrs_post_init__(self) -> None:
        assert self.suffix.startswith("."), (
            f"suffix must start with a dot: {self.suffix!r}"
        )

    @property
    def deps(self) -> list[Task]:
        return [self.source_directory_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source = fetch(self.source_directory_task.asset_id)
        assert isinstance(source, DirectoryAsset), (
            f"{self.source_directory_task.asset_id} is not a directory asset"
        )
        for path in sorted(source.path.iterdir()):
            if path.is_file() and path.suffix == self.suffix:
                shutil.copyfile(path, scratch_dir / path.name)
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(cls, asset_id: str, source_directory_task: Task, suffix: str) -> Self:
        source_meta = source_directory_task.meta
        assert isinstance(source_meta, ResultDirectoryMeta), (
            f"source_directory_task must carry ResultDirectoryMeta, got {type(source_meta).__name__}"
        )
        return cls(
            meta=GWASPlotDirectoryMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
            ),
            source_directory_task=source_directory_task,
            suffix=suffix,
        )
