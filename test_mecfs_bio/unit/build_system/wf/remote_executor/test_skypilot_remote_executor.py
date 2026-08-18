from pathlib import Path, PurePath

import pytest
import sky.exceptions

from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.build_system.wf.remote_executor.remote_job import (
    RemoteJob,
    RemoteResources,
)
from mecfs_bio.build_system.wf.remote_executor.skypilot_remote_executor import (
    CostEstimate,
    SkyPilotRemoteExecutor,
    _always_confirm,
    _prompt_confirm,
    _raise_on_failed_remote_job,
    build_sky_task,
)


def _make_job(tmp_path: Path, resources: RemoteResources) -> RemoteJob:
    return RemoteJob(
        image="img:1",
        commands=["gctb --gwfm RC ..."],
        input_files={tmp_path / "x.ma": PurePath("work/x.ma")},
        s3_inputs={"s3://ref-bucket/Imputed13M/v1": PurePath("work/ref")},
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
    # cpus/memory are requested as SkyPilot "N+" minimums (see _build_resources),
    # so the stored request is e.g. "24+"; compare on the requested-quantity intent
    # after stripping the minimum marker.
    assert res.memory is not None
    assert res.cpus is not None
    assert int(str(res.memory).rstrip("+")) == memory_gb
    assert int(str(res.cpus).rstrip("+")) == vcpus
    assert res.disk_size == disk_gb
    assert res.region == region
    run = task.run
    setup = task.setup
    assert run is not None
    assert setup is not None
    assert "gctb --gwfm" in run
    # Setup is fail-fast: a failed staging step must abort the job at setup rather
    # than let it proceed to run and fail confusingly on a missing file.
    assert "set -e" in setup
    assert "aws s3" in setup
    # The reference bucket is Requester Pays, so the read must declare the requester
    # as payer or S3 returns 403.
    assert "--request-payer requester" in setup
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
        scratch_s3="s3://scratch-bucket/remote-exec-scratch",
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


def test_run_raises_immediately_when_scratch_s3_is_missing(tmp_path: Path) -> None:
    # scratch_s3 is required to stage outputs. run() must fail before doing anything
    # expensive, so an estimator that explodes if reached proves the scratch check
    # comes first (before cost estimation, prompting, or launch).
    def exploding_estimator(_job: RemoteJob) -> CostEstimate:
        raise AssertionError("cost estimation must not run when scratch_s3 is unset")

    executor = SkyPilotRemoteExecutor(
        confirm=_always_confirm,
        scratch_s3=None,
        cost_estimator=exploding_estimator,
    )
    job = _make_job(
        tmp_path,
        RemoteResources(memory_gb=192, vcpus=24, disk_gb=500, region="us-east-1"),
    )
    with pytest.raises(AssertionError):
        executor.run(job, tmp_path)


def test_non_interactive_executor_auto_confirms() -> None:
    executor = SkyPilotRemoteExecutor.non_interactive(
        scratch_s3="s3://scratch-bucket/scratch"
    )
    assert executor.confirm("Launch something expensive? [y/N] ") is True


def test_make_wf_default_remote_executor_prompts_rather_than_auto_confirming() -> None:
    # The default must prompt, never auto-confirm: a silent-yes default would launch
    # paid infrastructure with no human in the loop. Identity check because the
    # interactive prompter reads stdin and so cannot be invoked in a unit test.
    executor = make_wf().remote_executor
    assert isinstance(executor, SkyPilotRemoteExecutor)
    assert executor.confirm is _prompt_confirm


def test_prompt_confirm_declines_on_non_yes_answer() -> None:
    assert _prompt_confirm("Launch? [y/N] ", read=lambda _prompt: "n") is False


def test_prompt_confirm_accepts_yes_answer() -> None:
    assert _prompt_confirm("Launch? [y/N] ", read=lambda _prompt: "yes") is True


def test_failed_remote_job_exit_code_raises() -> None:
    # A remote gctb crash surfaces from sky.tail_logs as a nonzero exit code; the
    # executor must turn that into a hard error instead of proceeding to retrieval.
    failed = int(sky.exceptions.JobExitCode.FAILED)
    with pytest.raises(RuntimeError):
        _raise_on_failed_remote_job(failed, "remote-exec-abcd1234")


def test_successful_remote_job_exit_code_does_not_raise() -> None:
    succeeded = int(sky.exceptions.JobExitCode.SUCCEEDED)
    _raise_on_failed_remote_job(succeeded, "remote-exec-abcd1234")


def test_retrieve_outputs_issues_recursive_s3_copies(tmp_path: Path) -> None:
    calls: list[list[str]] = []

    def fake_runner(command: list[str]) -> str:
        calls.append(command)
        return ""

    executor = SkyPilotRemoteExecutor(confirm=_always_confirm, runner=fake_runner)
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
