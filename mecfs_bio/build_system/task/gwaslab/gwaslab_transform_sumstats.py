"""
Contains a task for transforming gwaslab Sumstats objefts.
"""

from pathlib import Path

import gwaslab as gl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_sumstats import read_sumstats
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GwasLabTransformSpec,
    transform_gwaslab_sumstats,
)
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class GWASLabTransformSumstatsTask(Task):
    """
    Task to read a pickled GWASlab sumstats object, apply transformations to it, write it out as a new GWASLab sumstats object.
    """

    _meta: Meta
    source_sumstats_task: Task
    transform_spec: GwasLabTransformSpec

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def source_meta(self) -> Meta:
        return self.source_sumstats_task.meta

    @property
    def source_asset_id(self) -> AssetId:
        return self.source_meta.asset_id

    @property
    def deps(self) -> list["Task"]:
        return [self.source_sumstats_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        asset = fetch(self.source_asset_id)
        sumstats = read_sumstats(asset)
        sumstats = transform_gwaslab_sumstats(
            sumstats,
            spec=self.transform_spec,
        )
        out_path = scratch_dir / "pickled_sumstats.pickle"
        gl.dump_pickle(sumstats, path=out_path)
        return FileAsset(out_path)

    @classmethod
    def create_from_source_task(
        cls,
        source_tsk: Task,
        asset_id: str,
        spec: GwasLabTransformSpec,
    ):
        source_meta = source_tsk.meta
        assert isinstance(source_meta, GWASLabSumStatsMeta)
        meta = GWASLabSumStatsMeta(
            short_id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            sub_dir=source_meta.sub_dir,
        )
        return cls(
            meta=meta,
            source_sumstats_task=source_tsk,
            transform_spec=spec,
        )
