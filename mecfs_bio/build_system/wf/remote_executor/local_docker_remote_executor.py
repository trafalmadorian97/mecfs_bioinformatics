import shutil
import tempfile
from collections.abc import Callable
from pathlib import Path

from mecfs_bio.build_system.wf.remote_executor.base_remote_executor import (
    RemoteExecutor,
)
from mecfs_bio.build_system.wf.remote_executor.remote_job import RemoteJob
from mecfs_bio.util.subproc.run_command import execute_command


class LocalDockerRemoteExecutor(RemoteExecutor):
    """Runs a RemoteJob's container locally via Docker.

    Used as the system-test seam for the remote-executor abstraction: it
    exercises the same RemoteJob contract (input staging, container
    invocation, output retrieval) as a real cloud executor, but does all of
    it on the local machine with a local Docker daemon. RemoteResources is
    intentionally ignored, since there is no remote instance to size.
    """

    def __init__(self, runner: Callable[[list[str]], str] = execute_command) -> None:
        self._runner = runner

    def run(self, job: RemoteJob, local_output_dir: Path) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            host_work_dir: Path = Path(tmp_dir)
            for local_source, remote_dest in job.input_files.items():
                staged_dest: Path = host_work_dir / remote_dest
                staged_dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy(local_source, staged_dest)

            for s3_uri, remote_dest in job.s3_inputs.items():
                staged_dest = host_work_dir / remote_dest
                staged_dest.parent.mkdir(parents=True, exist_ok=True)
                self._runner(["aws", "s3", "cp", s3_uri, str(staged_dest)])

            inner_command: str = " && ".join(job.commands)
            docker_command: list[str] = [
                "docker",
                "run",
                "--rm",
                "-v",
                f"{host_work_dir}:/work",
                "-w",
                "/work",
                job.image,
                "bash",
                "-lc",
                f'"{inner_command}"',
            ]
            self._runner(docker_command)

            for output_file in job.output_files:
                staged_output: Path = host_work_dir / output_file
                assert staged_output.exists(), (
                    f"LocalDockerRemoteExecutor: expected container output "
                    f"{output_file} at {staged_output}, but it does not exist"
                )
                final_dest: Path = local_output_dir / output_file
                final_dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy(staged_output, final_dest)
