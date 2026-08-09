from abc import ABC, abstractmethod
from pathlib import Path

from mecfs_bio.build_system.wf.remote_executor.remote_job import RemoteJob


class RemoteExecutor(ABC):
    """Runs a RemoteJob's container commands somewhere and retrieves its outputs."""

    @abstractmethod
    def run(self, job: RemoteJob, local_output_dir: Path) -> None: ...
