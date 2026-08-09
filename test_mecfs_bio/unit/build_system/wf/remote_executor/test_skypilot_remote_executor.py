from pathlib import Path, PurePath

import pytest

from mecfs_bio.build_system.wf.remote_executor.remote_job import (
    RemoteJob,
    RemoteResources,
)
from mecfs_bio.build_system.wf.remote_executor.skypilot_remote_executor import (
    CostEstimate,
    SkyPilotRemoteExecutor,
    _prompt_confirm,
    build_sky_task,
)


def _make_job(tmp_path: Path, resources: RemoteResources) -> RemoteJob:
    return RemoteJob(
        image="img:1",
        commands=["gctb --gwfm RC ..."],
        input_files={tmp_path / "x.ma": PurePath("work/x.ma")},
        s3_inputs={"s3://b/Imputed13M/v1": PurePath("work/ref")},
        output_files=[PurePath("work/out")],
        resources=resources,
    )


def test_build_sky_task_reflects_resources_and_commands(tmp_path: Path) -> None:
    memory_gb = 192
    vcpus = 24
    disk_gb = 500
    region = "us-east-1"
    job = _make_job(
        tmp_path,
        RemoteResources(
            memory_gb=memory_gb, vcpus=vcpus, disk_gb=disk_gb, region=region
        ),
    )
    task = build_sky_task(job)
    res = list(task.resources)[0]
    # SkyPilot normalizes memory/cpus to strings (e.g. "192"); compare on the
    # requested-quantity intent rather than the exact attribute representation.
    assert res.memory is not None
    assert res.cpus is not None
    assert int(res.memory) == memory_gb
    assert int(res.cpus) == vcpus
    assert res.disk_size == disk_gb
    assert res.region == region
    run = task.run
    setup = task.setup
    assert run is not None
    assert setup is not None
    assert "gctb --gwfm" in run
    assert "aws s3" in setup
    # The container image is pulled in setup and run inside the container in run.
    assert "docker pull" in setup
    assert "docker run" in run


def test_run_prompt_shows_injected_cost_estimate_and_decline_aborts(
    tmp_path: Path,
) -> None:
    # The injected cost estimator stands in for the SkyPilot-optimizer call so the
    # unit test never contacts a cloud; declining the prompt must abort before any
    # launch, and the prompt must surface the concrete instance and rate.
    estimate = CostEstimate(
        usd_per_hour=1.53,
        instance_type="m6i.8xlarge",
        cloud="AWS",
        region="us-east-1",
    )
    seen_prompts: list[str] = []

    def decline(prompt: str) -> bool:
        seen_prompts.append(prompt)
        return False

    executor = SkyPilotRemoteExecutor(
        confirm=decline,
        cost_estimator=lambda _job: estimate,
    )
    job = _make_job(
        tmp_path,
        RemoteResources(memory_gb=192, vcpus=24, disk_gb=500, region="us-east-1"),
    )
    with pytest.raises(RuntimeError):
        executor.run(job, tmp_path)
    assert len(seen_prompts) == 1
    assert estimate.instance_type in seen_prompts[0]
    assert f"{estimate.usd_per_hour:.2f}" in seen_prompts[0]


def test_prompt_confirm_honours_assume_yes_env(monkeypatch) -> None:
    monkeypatch.setenv("REMOTE_EXEC_ASSUME_YES", "1")
    # The env override wins even if the injected reader would decline.
    assert _prompt_confirm("Launch? [y/N] ", read=lambda _prompt: "n") is True


def test_prompt_confirm_declines_on_non_yes_answer(monkeypatch) -> None:
    monkeypatch.delenv("REMOTE_EXEC_ASSUME_YES", raising=False)
    assert _prompt_confirm("Launch? [y/N] ", read=lambda _prompt: "n") is False


def test_prompt_confirm_accepts_yes_answer(monkeypatch) -> None:
    monkeypatch.delenv("REMOTE_EXEC_ASSUME_YES", raising=False)
    assert _prompt_confirm("Launch? [y/N] ", read=lambda _prompt: "yes") is True


def test_retrieve_outputs_issues_recursive_s3_copies(tmp_path: Path) -> None:
    calls: list[list[str]] = []

    def fake_runner(command: list[str]) -> str:
        calls.append(command)
        return ""

    executor = SkyPilotRemoteExecutor(runner=fake_runner)
    output_files = [PurePath("work/out"), PurePath("work/log.txt")]
    executor._retrieve_outputs(
        "s3://bucket/remote-exec-scratch/remote-exec-abcd1234",
        output_files,
        tmp_path,
    )

    assert len(calls) == len(output_files)
    for command, output_file in zip(calls, output_files, strict=True):
        assert command[:4] == ["aws", "s3", "cp", "--recursive"]
        source, dest = command[4], command[5]
        assert output_file.as_posix() in source
        assert str(tmp_path / output_file) in dest
    # Destination parent directories are created ahead of the copy.
    assert (tmp_path / "work").is_dir()
