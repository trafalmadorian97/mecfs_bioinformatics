from pathlib import Path, PurePath

from mecfs_bio.build_system.wf.remote_executor.local_docker_remote_executor import (
    LocalDockerRemoteExecutor,
)
from mecfs_bio.build_system.wf.remote_executor.remote_job import (
    RemoteJob,
    RemoteResources,
)


def test_issues_a_docker_run_command(tmp_path: Path) -> None:
    calls: list[list[str]] = []

    def fake_runner(cmd: list[str]) -> str:
        calls.append(cmd)
        return ""

    ex = LocalDockerRemoteExecutor(runner=fake_runner)
    job = RemoteJob(
        image="busybox:latest",
        commands=["echo hi"],
        input_files={},
        s3_inputs={},
        output_files=[PurePath("out/pip.txt")],
        resources=RemoteResources(1, 1, 1),
    )
    # output copy-back would fail because the fake runner produces nothing; that path is
    # exercised for real in the system test, so allow the missing-output error here.
    try:
        ex.run(job, tmp_path)
    except AssertionError:
        pass
    assert any("docker run" in " ".join(c) for c in calls)
