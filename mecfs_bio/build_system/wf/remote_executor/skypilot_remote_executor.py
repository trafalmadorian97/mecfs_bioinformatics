"""Production AWS remote executor built on the SkyPilot Python SDK.

This executor provisions a transient on-demand EC2 instance via SkyPilot, stages
the RemoteJob's inputs, runs its container command, retrieves the outputs, and
guarantees teardown of the instance (even on exception or Ctrl-C).

Assumptions and environment:

- AWS credentials are resolved through the standard AWS credential chain
  (environment variables, shared credentials file, or an instance profile);
  neither this module nor SkyPilot is passed explicit credentials. SkyPilot and
  the on-instance aws CLI both rely on that chain.
- On-demand instances only. Spot is intentionally not used (see design spec).
- Because a real launch runs for many hours and costs real money, run() asks for
  interactive confirmation before launching. Set the environment variable
  REMOTE_EXEC_ASSUME_YES=1 to skip the prompt in non-interactive contexts (CI
  smoke tests, batch runs).
- Output retrieval uses an S3 round-trip (see _retrieve_outputs): the on-instance
  run phase copies each output to a run-scoped S3 prefix, and run() downloads
  from that prefix. The scratch bucket/prefix is read from the environment
  variable REMOTE_EXEC_SCRATCH_S3 (e.g. "s3://my-bucket/remote-exec-scratch").
  The installed SkyPilot SDK exposes no clean, recursive, top-level file-download
  API (only log download and Storage-mount helpers), so the documented S3
  fallback is used.

The confirmation prompt shows a concrete on-demand price: the injected
cost_estimator asks SkyPilot's optimizer for the cheapest matching instance and
its hourly rate. The production estimator (estimate_cost_via_sky_optimize) needs
SkyPilot cloud access, so it is injected on construction and replaced with a fake
in unit tests; the pure helper build_sky_task is unit-tested directly, and the
live launch path is validated by a one-time manual maintainer smoke test.
"""

import os
import shlex
from collections.abc import Callable, Sequence
from pathlib import Path, PurePath
from uuid import uuid4

import sky
import structlog
from attrs import frozen

from mecfs_bio.build_system.wf.remote_executor.base_remote_executor import (
    RemoteExecutor,
)
from mecfs_bio.build_system.wf.remote_executor.remote_job import (
    RemoteJob,
    RemoteResources,
)
from mecfs_bio.util.subproc.run_command import execute_command

logger = structlog.get_logger()

_SECONDS_PER_HOUR = 3600

_SCRATCH_S3_ENV_VAR = "REMOTE_EXEC_SCRATCH_S3"
_ASSUME_YES_ENV_VAR = "REMOTE_EXEC_ASSUME_YES"


@frozen
class CostEstimate:
    """The concrete on-demand instance SkyPilot's optimizer would pick, and its rate.

    usd_per_hour is that instance's on-demand price for one hour; instance_type,
    cloud, and region identify the machine so the launch prompt can show what will
    actually be provisioned. Not a billing guarantee, but a catalog-derived figure
    rather than a hand-tuned model.
    """

    usd_per_hour: float
    instance_type: str
    cloud: str
    region: str | None


def _prompt_confirm(prompt: str, read: Callable[[str], str] = input) -> bool:
    """Return True if the user confirms the prompt (or REMOTE_EXEC_ASSUME_YES=1).

    The environment override lets non-interactive callers proceed without a TTY;
    otherwise only an answer beginning with "y" (case-insensitive) confirms. read
    is injected (defaulting to input) so tests can supply their own reader.
    """
    if os.environ.get(_ASSUME_YES_ENV_VAR) == "1":
        return True
    answer = read(prompt)
    return answer.strip().lower().startswith("y")


def _build_resources(resources: RemoteResources) -> sky.Resources:
    """Translate our RemoteResources into a SkyPilot Resources request on AWS."""
    return sky.Resources(
        cloud=sky.AWS(),
        cpus=resources.vcpus,
        memory=resources.memory_gb,
        disk_size=resources.disk_gb,
        region=resources.region,
    )


def estimate_cost_via_sky_optimize(job: RemoteJob) -> CostEstimate:
    """Ask SkyPilot's optimizer for the cheapest matching instance and its rate.

    Builds a resources-only probe task (no file mounts, so nothing is validated or
    synced), runs it through sky.optimize, and reads the chosen best_resources and
    its one-hour on-demand cost. This connects to the SkyPilot API server and needs
    cloud access, which is why it is injected into SkyPilotRemoteExecutor and
    replaced with a fake in unit tests.
    """
    probe = sky.Task(name="cost-estimate", run="true")
    probe.set_resources(_build_resources(job.resources))
    with sky.Dag() as dag:
        dag.add(probe)
    best = sky.get(sky.optimize(dag)).tasks[0].best_resources
    assert best is not None and best.instance_type is not None, (
        "sky.optimize returned no feasible instance for the requested resources"
    )
    return CostEstimate(
        usd_per_hour=best.get_cost(_SECONDS_PER_HOUR),
        instance_type=best.instance_type,
        cloud=str(best.cloud),
        region=best.region,
    )


def build_sky_task(job: RemoteJob, output_s3_prefix: str | None = None) -> sky.Task:
    """Build a sky.Task that runs the RemoteJob's container on a remote instance.

    input_files are mapped to SkyPilot file_mounts (remote destination -> local
    source). setup copies each s3_inputs prefix onto the instance and pulls the
    container image. run invokes the container over the mounted working
    directory. When output_s3_prefix is given, run also uploads each output to
    that run-scoped prefix so run() can download it afterwards.
    """
    file_mounts: dict[str, str] = {
        str(remote_dest): str(local_source)
        for local_source, remote_dest in job.input_files.items()
    }

    setup_lines: list[str] = []
    for s3_uri, remote_dest in job.s3_inputs.items():
        setup_lines.append(
            f"aws s3 cp --recursive {shlex.quote(s3_uri)} "
            f"{shlex.quote(str(remote_dest))}"
        )
    setup_lines.append(f"docker pull {shlex.quote(job.image)}")
    setup = "\n".join(setup_lines)

    inner_command = " && ".join(job.commands)
    run_parts: list[str] = [
        f"docker run --rm -v $(pwd):/work -w /work {shlex.quote(job.image)} "
        f"bash -lc {shlex.quote(inner_command)}"
    ]
    if output_s3_prefix is not None:
        prefix = output_s3_prefix.rstrip("/")
        for output_file in job.output_files:
            dest = f"{prefix}/{PurePath(output_file).as_posix()}"
            run_parts.append(
                f"aws s3 cp --recursive {shlex.quote(str(output_file))} "
                f"{shlex.quote(dest)}"
            )
    run = " && ".join(run_parts)

    task = sky.Task(name="remote-job", setup=setup, run=run)
    task.set_file_mounts(file_mounts)
    task.set_resources(_build_resources(job.resources))
    return task


class SkyPilotRemoteExecutor(RemoteExecutor):
    """Runs a RemoteJob on a transient AWS instance provisioned by SkyPilot.

    Confirms before launch (interactively, or via REMOTE_EXEC_ASSUME_YES=1),
    provisions
    an on-demand instance with an idle-autostop safety net, streams logs, pulls
    outputs back through an S3 scratch prefix, and always tears the cluster down.
    """

    def __init__(
        self,
        confirm: Callable[[str], bool] = _prompt_confirm,
        idle_minutes_to_autostop: int = 15,
        runner: Callable[[list[str]], str] = execute_command,
        cost_estimator: Callable[[RemoteJob], CostEstimate] = (
            estimate_cost_via_sky_optimize
        ),
    ) -> None:
        self._confirm = confirm
        self._idle_minutes_to_autostop = idle_minutes_to_autostop
        self._runner = runner
        self._cost_estimator = cost_estimator

    def run(self, job: RemoteJob, local_output_dir: Path) -> None:
        estimate = self._cost_estimator(job)
        prompt = (
            f"Launch {job.resources.vcpus} vCPU / {job.resources.memory_gb} GB on "
            f"{estimate.instance_type} ({estimate.cloud}/{estimate.region}) on-demand "
            f"(~${estimate.usd_per_hour:.2f}/hr)? [y/N] "
        )
        if not self._confirm(prompt):
            raise RuntimeError("Remote launch declined by user")

        scratch_root = os.environ.get(_SCRATCH_S3_ENV_VAR)
        assert scratch_root, (
            f"{_SCRATCH_S3_ENV_VAR} must be set to an s3:// scratch prefix so "
            f"remote outputs can be staged for download"
        )
        cluster = f"remote-exec-{uuid4().hex[:8]}"
        output_s3_prefix = f"{scratch_root.rstrip('/')}/{cluster}"
        task = build_sky_task(job, output_s3_prefix=output_s3_prefix)
        try:
            request_id = sky.launch(
                task,
                cluster_name=cluster,
                idle_minutes_to_autostop=self._idle_minutes_to_autostop,
                down=True,
            )
            job_id, _ = sky.get(request_id)
            sky.tail_logs(cluster, job_id, follow=True)
            self._retrieve_outputs(output_s3_prefix, job.output_files, local_output_dir)
        finally:
            # Guaranteed teardown even on exception / Ctrl-C. sky.down is async and
            # returns a request id; await it so a failed teardown surfaces instead of
            # silently leaving a paid instance running (idle-autostop is only a
            # backstop). Never let a teardown error mask the original exception.
            try:
                sky.get(sky.down(cluster))
            except Exception:
                logger.error(
                    "remote_teardown_failed",
                    cluster=cluster,
                    hint="Instance may still be running; tear it down manually "
                    "(e.g. `sky down <cluster>`) to stop incurring cost.",
                    exc_info=True,
                )

    def _retrieve_outputs(
        self,
        output_s3_prefix: str,
        output_files: Sequence[PurePath],
        local_output_dir: Path,
    ) -> None:
        """Download run outputs from the run-scoped S3 prefix into local_output_dir.

        SkyPilot's Python SDK offers no clean recursive file-download API, so the
        run phase already uploaded each output to output_s3_prefix; here we mirror
        that prefix back down with the aws CLI.
        """
        for output_file in output_files:
            source = (
                f"{output_s3_prefix.rstrip('/')}/{PurePath(output_file).as_posix()}"
            )
            dest = local_output_dir / output_file
            dest.parent.mkdir(parents=True, exist_ok=True)
            self._runner(
                [
                    "aws",
                    "s3",
                    "cp",
                    "--recursive",
                    shlex.quote(source),
                    shlex.quote(str(dest)),
                ]
            )
