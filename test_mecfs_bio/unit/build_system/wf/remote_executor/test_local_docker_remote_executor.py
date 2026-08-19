from pathlib import Path, PurePath

from mecfs_bio.build_system.wf.remote_executor.local_docker_remote_executor import (
    LocalDockerRemoteExecutor,
)
from mecfs_bio.build_system.wf.remote_executor.remote_job import (
    RemoteJob,
    RemoteResources,
)


def _mount_host_dir(docker_cmd: list[str]) -> Path:
    """Recover the host directory bind-mounted at /work from a docker run command."""
    spec = docker_cmd[docker_cmd.index("-v") + 1]
    return Path(spec.rsplit(":/work", 1)[0])


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


def test_s3_uri_source_is_fetched_via_aws(tmp_path: Path) -> None:
    """A non-local (s3://) source is staged through the injected runner as aws s3 cp."""
    calls: list[list[str]] = []

    def fake_runner(cmd: list[str]) -> str:
        calls.append(cmd)
        if "run" in cmd:  # simulate the container producing its declared output
            out = _mount_host_dir(cmd) / "work" / "out.txt"
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text("done")
        return ""

    ex = LocalDockerRemoteExecutor(runner=fake_runner)
    job = RemoteJob(
        image="busybox:latest",
        commands=["echo hi"],
        input_files={},
        s3_inputs={"s3://bucket/prefix/": PurePath("work/ref")},
        output_files=[PurePath("work/out.txt")],
        resources=RemoteResources(1, 1, 1),
    )
    ex.run(job, tmp_path)

    aws_calls = [c for c in calls if "aws" in c and "s3" in c]
    assert aws_calls, "expected an aws s3 cp call for the s3:// source"
    assert any("s3://bucket/prefix/" in c for c in aws_calls)


def test_local_dir_s3_input_staged_and_directory_output_retrieved(
    tmp_path: Path,
) -> None:
    """A local-directory s3_inputs source is copytree-staged (no aws), and a directory
    output is retrieved with its contents intact."""
    ref_src = tmp_path / "ref_src"
    (ref_src / "sub").mkdir(parents=True)
    (ref_src / "a.txt").write_text("alpha")
    (ref_src / "sub" / "b.txt").write_text("beta")

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    calls: list[list[str]] = []
    staged_ref_files: list[str] = []

    def fake_runner(cmd: list[str]) -> str:
        calls.append(cmd)
        if "run" in cmd:  # simulate the container: read staged inputs, write outputs
            host = _mount_host_dir(cmd)
            staged_ref = host / "work" / "ref"
            staged_ref_files.extend(
                sorted(
                    str(p.relative_to(staged_ref))
                    for p in staged_ref.rglob("*")
                    if p.is_file()
                )
            )
            out = host / "work" / "out"
            out.mkdir(parents=True, exist_ok=True)
            (out / "result.txt").write_text("gamma")
        return ""

    ex = LocalDockerRemoteExecutor(runner=fake_runner)
    job = RemoteJob(
        image="busybox:latest",
        commands=["echo hi"],
        input_files={},
        s3_inputs={str(ref_src): PurePath("work/ref")},
        output_files=[PurePath("work/out")],
        resources=RemoteResources(1, 1, 1),
    )
    ex.run(job, output_dir)

    # Local-directory source is staged locally, never through aws.
    assert not any("aws" in c for c in calls)
    # Its files (including nested) land under the staged dest.
    assert set(staged_ref_files) == {"a.txt", str(PurePath("sub") / "b.txt")}
    # The directory output is retrieved as a directory with contents intact.
    retrieved = output_dir / "work" / "out"
    assert retrieved.is_dir()
    assert (retrieved / "result.txt").read_text() == "gamma"
