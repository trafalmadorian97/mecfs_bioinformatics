"""
A task that fetches GWAS data from the Open Science data store
"""

import shlex
from collections.abc import Callable
from pathlib import Path

from attrs import frozen

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.util.download.verify import verify_hash
from mecfs_bio.util.subproc.run_command import execute_command_with_retries


@frozen
class OSFRetrievalTask(GeneratingTask):
    """
    A task that fetches GWAS data from the Open Science data store

    osfclient issues exactly one HTTP request per download and raises on any
    non-200 response, so a transient error from OSF's file storage aborts the
    whole fetch. The retries here supply the backoff osfclient lacks.

    The fetch is forced because osfclient creates the local file before it
    starts downloading into it: after a failed attempt an empty file is left
    behind, and without --force osfclient refuses to overwrite it, which would
    make every subsequent retry fail for a second, unrelated reason.
    """

    meta: GWASSummaryDataFileMeta
    osf_project_id: str
    md5_hash: str | None = None
    run_command: Callable[[list[str]], str] = execute_command_with_retries

    def __attrs_post_init__(self):
        # project_path is optional on GWASSummaryDataFileMeta, since a GWAS that
        # does not come from OSF has no path within an OSF project. This task
        # cannot do anything without one, so reject it at construction rather
        # than at execution, which may be a long build later.
        assert self.meta.project_path is not None, (
            f"{type(self).__name__} needs a project_path to fetch from OSF project "
            f"{self.osf_project_id}, but meta {self.meta.asset_id} has none"
        )

    @property
    def deps(self) -> list[Task]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> FileAsset:
        tmp_dst = scratch_dir / "tmp"
        self.run_command(
            [
                "pixi",
                "r",
                "osf",
                "-p",
                self.osf_project_id,
                "fetch",
                "--force",
                shlex.quote(str(self.meta.project_path)),
                shlex.quote(str(tmp_dst)),
            ]
        )
        verify_hash(tmp_dst, self.md5_hash)

        return FileAsset(
            tmp_dst,
        )
