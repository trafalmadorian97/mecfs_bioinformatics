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
- Because a real launch runs for many hours and costs real money, run() asks the
  injected confirm callable before launching. Construct via SkyPilotRemoteExecutor
  .interactive() to prompt a human (the default), or .non_interactive() to
  auto-approve in a pipeline with no user present. The confirm callable has no
  default on the class, so the choice must be made explicitly at construction.
- Output retrieval uses an S3 round-trip (see _retrieve_outputs): the on-instance
  run phase copies each output to a run-scoped S3 prefix, and run() downloads
  from that prefix. The scratch bucket/prefix is the executor's scratch_s3
  attribute (e.g. "s3://my-bucket/remote-exec-scratch"), supplied on construction
  (in production, from default_runner_config.yaml). The installed SkyPilot SDK
  exposes no clean, recursive, top-level file-download API (only log download and
  Storage-mount helpers), so the documented S3 fallback is used.

The confirmation prompt shows a concrete on-demand price: the injected
cost_estimator asks SkyPilot's optimizer for the cheapest matching instance and
its hourly rate. The production estimator (estimate_cost_via_sky_optimize) needs
SkyPilot cloud access, so it is injected on construction and replaced with a fake
in unit tests; the pure helper build_sky_task is unit-tested directly, and the
live launch path is validated by a one-time manual maintainer smoke test.
"""

import shlex
from collections.abc import Callable, Sequence
from pathlib import Path, PurePath
from uuid import uuid4

import sky
import sky.exceptions
import structlog
from attrs import frozen

from mecfs_bio.build_system.wf.remote_executor.base_remote_executor import (
    RemoteExecutor,
)
from mecfs_bio.build_system.wf.remote_executor.remote_job import (
    RemoteJob,
    RemoteResources,
)
from mecfs_bio.util.format_verify.s3_uri import is_valid_s3_uri
from mecfs_bio.util.subproc.run_command import execute_command

logger = structlog.get_logger()

_SECONDS_PER_HOUR = 3600


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
    """Return True if the user answers the prompt affirmatively.

    Only an answer beginning with "y" (case-insensitive) confirms. read is injected
    (defaulting to input) so tests can supply their own reader. This is the confirm
    callable wired by SkyPilotRemoteExecutor.interactive().
    """
    answer = read(prompt)
    return answer.strip().lower().startswith("y")


def _always_confirm(prompt: str) -> bool:
    """Confirm callable that approves unconditionally, for non-interactive pipelines.

    Wired by SkyPilotRemoteExecutor.non_interactive() so a batch/CI run with no user
    present proceeds without prompting. Because this auto-approves a real, paid
    launch, it must be selected explicitly (never a default).
    """
    return True


def _raise_on_failed_remote_job(exit_code: int, cluster: str) -> None:
    """Raise if a finished remote job did not succeed.

    sky.tail_logs(follow=True) returns the job's exit code (JobExitCode.SUCCEEDED
    is 0; any other value means the job failed or was cancelled) rather than raising
    on failure. Callers must not swallow that: a nonzero code here means the gctb
    container crashed on the remote instance, and we surface it as a hard error so a
    remote failure never masquerades as a successful (but output-less) run.
    """
    if exit_code != int(sky.exceptions.JobExitCode.SUCCEEDED):
        raise RuntimeError(
            f"Remote job on cluster {cluster} failed with exit code {exit_code} "
            f"(sky.exceptions.JobExitCode); see the streamed logs above for the "
            f"container error."
        )


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
        # This executor provisions AWS and stages inputs with `aws s3 cp`, so every
        # declared input must be an s3:// URI. Assert it here rather than let a foreign
        # scheme (e.g. a gs:// URI from a different ObjectStore) fail deep inside the
        # remote aws CLI: the input store's cloud must match the executor's cloud, and
        # that constraint is otherwise unexpressed in code.
        assert is_valid_s3_uri(s3_uri), (
            f"SkyPilotRemoteExecutor provisions AWS and can only stage valid s3:// "
            f"inputs; got {s3_uri!r}. The reference ObjectStore's cloud must match "
            f"the remote executor's cloud."
        )
        # --request-payer requester is required to read the Requester Pays reference
        # bucket: the downloading instance (not the bucket owner) is billed for the
        # transfer, and the request is rejected with 403 without this flag.
        setup_lines.append(
            f"aws s3 cp --recursive --request-payer requester "
            f"{shlex.quote(s3_uri)} {shlex.quote(str(remote_dest))}"
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


@frozen
class SkyPilotRemoteExecutor(RemoteExecutor):
    """Runs a RemoteJob on a transient AWS instance provisioned by SkyPilot.

    Confirms before launch via the injected confirm callable, provisions an
    on-demand instance with an idle-autostop safety net, streams logs, pulls
    outputs back through an S3 scratch prefix, and always tears the cluster down.

    confirm has no class default: construct through the interactive() classmethod to
    prompt a human, or non_interactive() to auto-approve in a pipeline with no user
    present. Forcing that choice keeps a paid, unattended launch from ever being the
    silent default.

    scratch_s3 is the s3:// prefix under which remote outputs are staged for download
    (see run and _retrieve_outputs). It has no sensible default (every user supplies
    their own bucket), so it defaults to None and run() fails fast if it is still
    None at launch time.

    Frozen: the executor holds only its injected collaborators (confirm, runner,
    cost_estimator) and configuration (scratch_s3, autostop); none is reassigned
    after construction.
    """

    confirm: Callable[[str], bool]
    scratch_s3: str | None = None
    idle_minutes_to_autostop: int = 15
    runner: Callable[[list[str]], str] = execute_command
    cost_estimator: Callable[[RemoteJob], CostEstimate] = estimate_cost_via_sky_optimize

    def __attrs_post_init__(self) -> None:
        # scratch_s3 may be None (a runner without remote-exec configured); run()
        # fails fast at launch in that case. But if one is supplied, catch a
        # malformed prefix now rather than deep inside the on-instance aws CLI.
        assert self.scratch_s3 is None or is_valid_s3_uri(self.scratch_s3), (
            f"scratch_s3 {self.scratch_s3!r} is not a valid s3 URI"
        )

    @classmethod
    def interactive(cls, scratch_s3: str | None = None) -> "SkyPilotRemoteExecutor":
        """Build an executor that prompts a human to confirm each launch."""
        return cls(confirm=_prompt_confirm, scratch_s3=scratch_s3)

    @classmethod
    def non_interactive(cls, scratch_s3: str | None = None) -> "SkyPilotRemoteExecutor":
        """Build an executor that auto-approves launches, for unattended pipelines."""
        return cls(confirm=_always_confirm, scratch_s3=scratch_s3)

    def run(self, job: RemoteJob, local_output_dir: Path) -> None:
        # Fail fast, before cost estimation or any launch: staging outputs is
        # impossible without a scratch prefix, so there is no point provisioning
        # (or prompting to pay for) an instance whose outputs cannot be retrieved.
        assert self.scratch_s3 is not None, (
            "SkyPilotRemoteExecutor.scratch_s3 is not set; supply an s3:// scratch "
            "prefix (in production via default_runner_config.yaml) so remote outputs "
            "can be staged for download."
        )
        estimate = self.cost_estimator(job)
        prompt = (
            f"Launch {job.resources.vcpus} vCPU / {job.resources.memory_gb} GB on "
            f"{estimate.instance_type} ({estimate.cloud}/{estimate.region}) on-demand "
            f"(~${estimate.usd_per_hour:.2f}/hr)? [y/N] "
        )
        if not self.confirm(prompt):
            raise RuntimeError("Remote launch declined by user")

        cluster = f"remote-exec-{uuid4().hex[:8]}"
        output_s3_prefix = f"{self.scratch_s3.rstrip('/')}/{cluster}"
        task = build_sky_task(job, output_s3_prefix=output_s3_prefix)
        try:
            request_id = sky.launch(
                task,
                cluster_name=cluster,
                idle_minutes_to_autostop=self.idle_minutes_to_autostop,
                down=True,
            )
            job_id, _ = sky.get(request_id)
            exit_code = sky.tail_logs(cluster, job_id, follow=True)
            assert isinstance(exit_code, int), (
                "sky.tail_logs(follow=True) must return the job exit code"
            )
            _raise_on_failed_remote_job(exit_code, cluster)
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
            self.runner(
                [
                    "aws",
                    "s3",
                    "cp",
                    "--recursive",
                    shlex.quote(source),
                    shlex.quote(str(dest)),
                ]
            )
