from pathlib import Path
from typing import Sequence

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    PhenotypeInfo,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class LavaPhenotypeDataSource:
    """
    A source Task providing tabular Gwas summary statistics pertaining to a phenotype
    columns are assumed to be in GWASLAB format
    (see: mecfs_bio/constants/gwaslab_constants.py)
    will be copied and converted to format expected by lava

    """

    task: Task
    alias: str
    pipe: DataProcessingPipe = IdentityPipe()
    sample_info: PhenotypeInfo | None = None

    @property
    def asset_id(self) -> AssetId:
        return self.task.asset_id


@frozen
class LDReferenceInfo:
    ld_ref_task: Task
    filename_prefix: str


@frozen
class LavaTask(Task):
    """
    Given a locus definition file does the following for each locus in the file:
    - estimates heritability all phenotypes  using LAVA
    - For each pair of phenotypes, if the heritabilities are both sufficiently significant, calculates genetic correlation at the locus using LAVA


    """

    _meta: Meta
    sources: Sequence[LavaPhenotypeDataSource]
    ld_reference_info: LDReferenceInfo
    lava_locus_definitions_task: Task
    ct_ldsc_task_for_overlap: Task

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [item.task for item in self.sources] + [
            self.ld_reference_info.ld_ref_task
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        # todo
        return DirectoryAsset(scratch_dir)
