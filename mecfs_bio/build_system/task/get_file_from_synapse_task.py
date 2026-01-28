from pathlib import Path

import structlog
import synapseclient
import synapseutils
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.getLogger()


@frozen
class GetFileFromSynapseTask(Task):
    """
    Task to get a file from synapse.org (see: https://docs.synapse.org/synapse-docs/faq)
    May require authentication.  See Getting Started.
    """

    _meta: Meta
    synid: str
    expected_filename: str

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        syn = synapseclient.login()
        files = synapseutils.syncFromSynapse(syn, self.synid, path=str(scratch_dir))
        logger.debug(f"downloaded: {files}")
        assert len(files) == 1
        assert files[0].name == self.expected_filename
        return FileAsset(Path(files[0].path))
