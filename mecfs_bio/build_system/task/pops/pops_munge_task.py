"""
Task to munge a directory of raw gene features into the chunked matrix format POPs
consumes.

This wraps POPs' munge_feature_directory.py. It produces a directory containing
<prefix>.mat.<i>.npy, <prefix>.cols.<i>.txt, and <prefix>.rows.txt, which is then
passed to PopsRunTask via --feature_mat_prefix.
"""

from pathlib import Path

import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pops.pops_utils import (
    GENE_ANNOT_RELATIVE_PATH,
    MUNGE_SCRIPT_NAME,
    MUNGED_FEATURES_PREFIX_NAME,
    NanPolicy,
    invoke_pops_script,
)
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()


@frozen
class PopsMungeFeatureDirectoryTask(Task):
    """Munge raw POPs gene features into chunked matrices.

    raw_features_task may produce either a DirectoryAsset (a directory of feature
    files, e.g. the example features_raw directory) or a FileAsset (a single feature
    table, e.g. the gunzipped full feature collection); in the latter case the file
    is staged into a directory since munge_feature_directory.py processes every file
    in --feature_dir.
    """

    meta: Meta
    pops_source_task: Task
    raw_features_task: Task
    max_cols: int = 500
    nan_policy: NanPolicy = "raise"
    save_prefix_name: str = MUNGED_FEATURES_PREFIX_NAME

    @property
    def deps(self) -> list["Task"]:
        return [self.pops_source_task, self.raw_features_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self.pops_source_task.asset_id)
        features_asset = fetch(self.raw_features_task.asset_id)
        assert isinstance(source_asset, DirectoryAsset)
        source_dir = source_asset.path

        feature_dir = self._resolve_feature_dir(features_asset, scratch_dir)
        gene_annot_path = source_dir / GENE_ANNOT_RELATIVE_PATH
        out_dir = scratch_dir / "munged_features"
        out_dir.mkdir(parents=True, exist_ok=True)
        save_prefix = out_dir / self.save_prefix_name

        args = [
            "--gene_annot_path",
            str(gene_annot_path),
            "--feature_dir",
            str(feature_dir),
            "--save_prefix",
            str(save_prefix),
            "--max_cols",
            str(self.max_cols),
            "--nan_policy",
            self.nan_policy,
        ]
        invoke_pops_script(source_dir, MUNGE_SCRIPT_NAME, args)
        return DirectoryAsset(out_dir)

    def _resolve_feature_dir(self, features_asset: Asset, scratch_dir: Path) -> Path:
        if isinstance(features_asset, DirectoryAsset):
            return features_asset.path
        assert isinstance(features_asset, FileAsset)
        # munge_feature_directory.py reads every file in --feature_dir, so stage the
        # single feature file into a dedicated directory. Symlink to avoid copying a
        # potentially multi-GB feature collection.
        feature_dir = scratch_dir / "feature_dir"
        feature_dir.mkdir(parents=True, exist_ok=True)
        staged = feature_dir / features_asset.path.name
        staged.symlink_to(features_asset.path)
        return feature_dir

    @classmethod
    def create(
        cls,
        asset_id: str,
        pops_source_task: Task,
        raw_features_task: Task,
        max_cols: int = 500,
        nan_policy: NanPolicy = "raise",
    ) -> "PopsMungeFeatureDirectoryTask":
        return cls(
            meta=SimpleDirectoryMeta(id=AssetId(asset_id)),
            pops_source_task=pops_source_task,
            raw_features_task=raw_features_task,
            max_cols=max_cols,
            nan_policy=nan_policy,
        )
