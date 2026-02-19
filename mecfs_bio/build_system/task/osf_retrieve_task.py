"""
A task that fetches GWAS data from the Open Science data store
"""

import shlex
from pathlib import Path

import invoke
from attrs import frozen

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.util.download.verify import verify_hash


@frozen
class OSFRetrievalTask(GeneratingTask):
    """
    A task that fetches GWAS data from the Open Science data store
    """

    _meta: GWASSummaryDataFileMeta
    osf_project_id: str
    md5_hash: str | None = None

    @property
    def meta(self) -> GWASSummaryDataFileMeta:
        return self._meta

    @property
    def deps(self) -> list[Task]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> FileAsset:
        tmp_dst = scratch_dir / "tmp"

        @invoke.task
        def fetch_osf(c):
            c.run(
                f"pixi r osf -p {self.osf_project_id} fetch {str(shlex.quote(str(self._meta.project_path)))} {str(tmp_dst)}"
            )

        fetch_osf(invoke.Context())
        verify_hash(tmp_dst, self.md5_hash)

        return FileAsset(
            tmp_dst,
        )
