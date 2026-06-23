from pathlib import Path

import gwaslab as gl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import GWASLabSumStatsMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_sumstats import read_sumstats
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class GWASLabRSIDAssignmentOptions:
    threads:int=6


@frozen
class GWASLabAssignRSIDSViaSweepTask(Task):
    """
    Task to use GWASLab's sweep functionality to assign RSIDs
    """
    sumstats_task: Task
    vcf_dir_task: Task
    vcf_filename:str
    meta: Meta
    options: GWASLabRSIDAssignmentOptions



    @property
    def source_sumstats_asset_id(self) -> AssetId:
        return self.sumstats_task.asset_id



    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        asset = fetch(self.source_sumstats_asset_id)
        sumstats = read_sumstats(asset)
        vcf_asset = fetch(self.vcf_dir_task.asset_id)
        assert isinstance(vcf_asset, DirectoryAsset)
        vcf_path = vcf_asset.path/ self.vcf_filename
        sumstats.assign_rsid2(
            vcf_path=str(vcf_path),
            thead=self.options.threads
        )


    @property
    def deps(self) -> list["Task"]:
        return [self.sumstats_task, self.vcf_dir_task]


    @classmethod
    def create(cls,
        id: str,
        sumstats_task: Task,
        vcf_dir_task: Task,
        vcf_filename: str,
        options: GWASLabRSIDAssignmentOptions = GWASLabRSIDAssignmentOptions()
    ):
        source_meta=sumstats_task.meta
        assert isinstance(source_meta, GWASLabSumStatsMeta)
        meta = GWASLabSumStatsMeta(
           id=AssetId(id),
            trait=source_meta.trait,
            project=source_meta.project,
            sub_dir=source_meta.sub_dir,
        )
        return cls(
            sumstats_task=sumstats_task,
            vcf_dir_task=vcf_dir_task,
            vcf_filename=vcf_filename,
            meta=meta,
            options=options
        )
