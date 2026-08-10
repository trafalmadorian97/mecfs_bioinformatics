from collections.abc import Mapping
from pathlib import Path, PurePath

from mecfs_bio.build_system.wf.remote_executor.base_remote_executor import (
    RemoteExecutor,
)
from mecfs_bio.build_system.wf.remote_executor.remote_job import RemoteJob


class FakeRemoteExecutor(RemoteExecutor):
    def __init__(self, stub_outputs: Mapping[PurePath, str] | None = None) -> None:
        self._stub_outputs: dict[PurePath, str] = dict(stub_outputs or {})
        self.last_job: RemoteJob | None = None

    def run(self, job: RemoteJob, local_output_dir: Path) -> None:
        self.last_job = job
        # Declared output_files entries are directories (the real executors retrieve
        # each via aws s3 cp --recursive), so create them as directories.
        for out in job.output_files:
            (local_output_dir / out).mkdir(parents=True, exist_ok=True)
        # Materialize any stub files within local_output_dir.
        for rel_path, content in self._stub_outputs.items():
            dest = local_output_dir / rel_path
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_text(content)
