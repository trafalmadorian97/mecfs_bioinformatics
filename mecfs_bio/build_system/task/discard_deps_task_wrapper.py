"""
Save disk space materializing dependencies of the wrapped task in a temporary directory.
"""

import tempfile
from pathlib import Path

import attrs
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class DiscardDepsWrapper(Task):
    """
    Save disk space materializing dependencies of the wrapped task in a temporary directory.
    """

    _inner: Task

    @property
    def meta(self) -> Meta:
        return attrs.evolve(
            self._inner.meta, id=AssetId(self._inner.meta.id + "_wrapped_remove_deps")
        )

    @property
    def deps(self) -> list["Task"]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        out_path = scratch_dir / "out"
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            info_store = temp_path / "info_store.yaml"

            asset_root = temp_path / "asset_store"
            runner = SimpleRunner(
                tracer=ImoHasher.with_xxhash_128(),
                info_store=info_store,
                asset_root=asset_root,
            )
            result = runner.run([self._inner])
            result_asset = result[self._inner.asset_id]
            if isinstance(result_asset, FileAsset):
                is_file = True
            else:
                is_file = False
                assert isinstance(result_asset, DirectoryAsset)
            result_path = result_asset.path
            result_path.rename(out_path)
        if is_file:
            return FileAsset(out_path)
        return DirectoryAsset(out_path)
