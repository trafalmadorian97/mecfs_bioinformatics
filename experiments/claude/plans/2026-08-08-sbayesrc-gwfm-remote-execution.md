# SBayesRC GWFM Remote Execution — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run SBayesRC genome-wide fine-mapping (GWFM) for a trait on a transient AWS instance while the build graph runs locally, with results flowing back as a normal cached build asset.

**Architecture:** A new `RemoteExecutor` capability on the build system's `WF` object dispatches a self-contained container command (the `gctb` GWFM run) to a SkyPilot-provisioned instance; small per-run inputs are shipped from the laptop, the large shared reference bundle is pulled from S3 on the instance, and outputs are retrieved into the local scratch dir so the verifying-trace rebuilder caches them normally. The fine-mapping Task only builds a `RemoteJob` — it is agnostic to the provisioning mechanism.

**Tech Stack:** Python 3.12, attrs (`@frozen`), polars, boto3 (S3), SkyPilot (AWS), Docker, pytest, pixi task runner, structlog.

## Global Constraints

- All commands run via `pixi r <command>`; Python via `pixi r python <script>`. (CLAUDE.md)
- **Type-annotate every function/method/test parameter and return** (including `tmp_path: Path`, fake-runner `cmd: list[str]`, `run(self, job: RemoteJob, local_output_dir: Path) -> None`). No bare parameters.
- **After each task**, run `pixi r invoke green 2>&1 | tee experiments/claude/logs/green_task<N>.log`, then open the logfile and confirm it reports a pass before committing — tailing shows only R console noise, and testmon can print "no tests ran", so exit 0 alone is not proof. (project memory)
- **Do not edit anything under `docs/`** — that tree is the human-written published mkdocs site. Put any plan/design/setup notes under `experiments/claude/`.
- In docstrings: no backticks around inline code, no RST. (CLAUDE.md)
- Prefer `Path` for filesystem paths, `PurePath` for base-relative paths; convert to `str` only at serialization boundaries. (project memory)
- Prefer polars over pandas; return `@frozen` objects with named attributes, not bare tuples. (project memory)
- Never monkeypatch/mock — inject dependencies (executor, object store, subprocess runner) as params with production defaults; tests pass their own. Shell out via `execute_command` (`mecfs_bio/util/subproc/run_command.py`), not raw `subprocess`. (project memory)
- Enforce cross-field invariants in `__attrs_post_init__` (fail fast). Do NOT then unit-test those trivial guards or other constant/field-presence facts — avoid brittle/cargo-cult tests. (project memory)
- Type small fixed-value string params as `Literal` aliases; test at the Task level where practical. (project memory)
- On-demand instances only. No `use_spot` in any interface unless the recon finds explicit GCTB MCMC-resumption evidence.
- Reference-bundle S3 objects: One Zone-IA storage class, bucket region == compute region.
- GCTB is MIT-licensed; the published image must include GCTB's LICENSE/copyright notice.

---

## File Structure

New package `mecfs_bio/build_system/wf/remote_executor/`:
- `remote_job.py` — `RemoteResources`, `RemoteJob` frozen types.
- `base_remote_executor.py` — `RemoteExecutor` ABC.
- `fake_remote_executor.py` — `FakeRemoteExecutor` (records job; writes stub outputs) for unit tests.
- `local_docker_remote_executor.py` — `LocalDockerRemoteExecutor` (runs the identical container locally) for the system test.
- `skypilot_remote_executor.py` — `SkyPilotRemoteExecutor` (production AWS) + pure helpers.

New object-store capability `mecfs_bio/build_system/wf/object_store/`:
- `base_object_store.py` — `ObjectStore` ABC + `ObjectHead` frozen result.
- `s3_object_store.py` — `S3ObjectStore` (boto3, One Zone-IA).
- `fake_object_store.py` — `FakeObjectStore` (in-memory) for unit tests.

GWFM domain (`mecfs_bio/build_system/task/sbayesrc/`):
- `gctb_gwfm_constants.py` — pinned versions, URLs, bundle manifest, CLI templates, resource/disk sizing, marker JSON keys. (Only values the code consumes; exploration findings live in `experiments/claude/`.)
- `stage_gwfm_reference_task.py` — `StageGwfmReferenceTask`.
- `sumstats_to_cojo_ma_task.py` — `SumstatsToCojoMaTask`.
- `gctb_fine_map_task.py` — `GctbFineMapTask`.

Wiring / config / infra:
- `mecfs_bio/build_system/wf/base_wf.py` — add `remote_executor`, `object_store` to `WF` + `make_wf`.
- `mecfs_bio/build_system/runner/simple_runner.py:115` — pass the two capabilities into `make_wf`.
- `mecfs_bio/analysis/runner/default_runner.py` — read remote-exec config (region, bucket, image); lazily construct the SkyPilot executor.
- `docker/gctb/Dockerfile`, `docker/gctb/LICENSE-GCTB` — minimal GCTB image.
- `tasks.py` — `build_push_gctb_image` invoke task.

Tests mirror under `test_mecfs_bio/unit/build_system/...` and one `test_mecfs_bio/system/test_gctb_gwfm.py`.

---

## Task 1: Reconnaissance — pin constants + record findings

No unit tests (testing plain constants is bloat). Deliverable: a constants module holding only what code consumes, plus a findings doc for the exploratory facts.

**Files:**
- Create: `mecfs_bio/build_system/task/sbayesrc/gctb_gwfm_constants.py`
- Create: `experiments/claude/design_specs/2026-08-08-gwfm-recon-findings.md`

**Interfaces:**
- Produces constants consumed by later tasks: `GCTB_VERSION: str`, `GCTB_BINARY_URL: str`, `GCTB_BINARY_SHA256: str`, `GWFM_REFERENCE_VERSION: str`, `ReferenceBundleFile` (frozen: `filename: str`, `source_url: str`, `size_bytes: int`, `sha256: str | None`), `GWFM_REFERENCE_BUNDLE: tuple[ReferenceBundleFile, ...]`, CLI templates `GCTB_MAKE_LDM_EIGEN_TEMPLATE`/`GCTB_GWFM_TEMPLATE`/`GCTB_CS_TEMPLATE: str`, `DEFAULT_MEMORY_GB`/`DEFAULT_VCPUS`/`DEFAULT_DISK_GB: int`, and marker JSON keys `MARKER_VERSION_KEY`/`MARKER_PREFIX_KEY`/`MARKER_FILES_KEY: str`.

- [ ] **Step 1: Investigate; write findings to `experiments/claude/design_specs/2026-08-08-gwfm-recon-findings.md`** (one source link each):
  1. gctbhub GWFM reference files + directory URLs under `.../resources/GWFM/LD/Imputed13M/`, and the annotation + gene-map location (check `.../resources/v2.0/Annotation/` and the GWFM tutorial). Record each file's `Content-Length` (HEAD).
  2. Binary URL `https://gctbhub.cloud.edu.au/software/gctb/download/gctb_2.5.5_Linux.zip`; download once, `sha256sum`, record.
  3. Exact GWFM CLI (the three `gctb` calls) and **whether `--gwfm RC` consumes the precomputed `eigen/` directly** (this decides whether Task 10 emits `--make-ldm-eigen` — a one-time authoring decision, recorded here, NOT a shipped boolean).
  4. Whether `gctb` GWFM supports checkpoint/resume (search repo/docs for resume/`.rds`/checkpoint). Record; if it does, spot may be revisited later — otherwise it stays out of scope.
  5. VM disk sizing (reference zip + unzip + eigen intermediates) → the value for `DEFAULT_DISK_GB` (~500).

- [ ] **Step 2: Write `gctb_gwfm_constants.py`** with the code-consumed values only. Bundle `sha256` fields start `None` (filled + committed after first real staging in Task 8). Shape:

```python
from attrs import frozen

GCTB_VERSION: str = "2.5.5"
GCTB_BINARY_URL: str = "https://gctbhub.cloud.edu.au/software/gctb/download/gctb_2.5.5_Linux.zip"
GCTB_BINARY_SHA256: str = "<recorded in Step 1.2>"

GWFM_REFERENCE_VERSION: str = "Imputed13M/v1"


@frozen
class ReferenceBundleFile:
    filename: str
    source_url: str
    size_bytes: int
    sha256: str | None  # None until first staging records it, then committed


GWFM_REFERENCE_BUNDLE: tuple[ReferenceBundleFile, ...] = (
    ReferenceBundleFile("ukbEUR_13M_FullLDM.zip", "https://.../ukbEUR_13M_FullLDM.zip", 0, None),
    # eigen files, annotation, gene-map ...
)

# Fields filled at run time: {ldm} {ma} {annot} {gene_map} {out} {threads} {pip} {pep} {pwld}
GCTB_MAKE_LDM_EIGEN_TEMPLATE: str = "gctb --ldm {ldm} --gwas-summary {ma} --make-ldm-eigen --thread {threads} --out {out}"
GCTB_GWFM_TEMPLATE: str = "gctb --gwfm RC --ldm-eigen {ldm} --gwas-summary {ma} --annot {annot} --gene-map {gene_map} --thread {threads} --out {out}"
GCTB_CS_TEMPLATE: str = "gctb --cs --pwld-file {pwld} --pip {pip} --pep {pep} --gene-map {gene_map} --mcmc-samples {out} --out {out}"

DEFAULT_MEMORY_GB: int = 192
DEFAULT_VCPUS: int = 24
DEFAULT_DISK_GB: int = 500

MARKER_VERSION_KEY: str = "version"
MARKER_PREFIX_KEY: str = "s3_prefix"
MARKER_FILES_KEY: str = "files"
```

- [ ] **Step 3: Run green** → `pixi r invoke green 2>&1 | tee experiments/claude/logs/green_task1.log`; confirm pass in the log.
- [ ] **Step 4: Commit.** `git add mecfs_bio/build_system/task/sbayesrc/gctb_gwfm_constants.py experiments/claude/design_specs/2026-08-08-gwfm-recon-findings.md && git commit -m "feat(gwfm): pin GCTB GWFM constants + record recon findings"`

---

## Task 2: `RemoteResources` and `RemoteJob` types

Pure type definitions with fail-fast invariants; no dedicated unit test (would be cargo-cult). Consumers exercise them.

**Files:**
- Create: `mecfs_bio/build_system/wf/remote_executor/remote_job.py`

**Interfaces:**
- Produces: `RemoteResources(memory_gb: int, vcpus: int, disk_gb: int, region: str | None = None)`; `RemoteJob(image: str, commands: Sequence[str], input_files: Mapping[Path, PurePath], s3_inputs: Mapping[str, PurePath], output_files: Sequence[PurePath], resources: RemoteResources)`.

- [ ] **Step 1: Implement.**

```python
from collections.abc import Mapping, Sequence
from pathlib import Path, PurePath

from attrs import frozen

# No `use_spot` / `instance_type` (see design spec). If a named-instance override is
# added later, refactor RemoteJob.resources to `RemoteResources | ExplicitInstance`
# rather than adding optional fields, keeping invalid states unrepresentable.


@frozen
class RemoteResources:
    memory_gb: int
    vcpus: int
    disk_gb: int
    region: str | None = None

    def __attrs_post_init__(self) -> None:
        assert self.memory_gb > 0 and self.vcpus > 0 and self.disk_gb > 0, (
            "RemoteResources fields must be positive"
        )


@frozen
class RemoteJob:
    image: str
    commands: Sequence[str]
    input_files: Mapping[Path, PurePath]
    s3_inputs: Mapping[str, PurePath]
    output_files: Sequence[PurePath]
    resources: RemoteResources

    def __attrs_post_init__(self) -> None:
        assert self.commands, "RemoteJob.commands must be non-empty"
        assert self.output_files, "RemoteJob.output_files must be non-empty"
```

- [ ] **Step 2: Run green** → `... | tee experiments/claude/logs/green_task2.log`; confirm pass.
- [ ] **Step 3: Commit.** `git commit -m "feat(remote): RemoteJob/RemoteResources types"`

---

## Task 3: `RemoteExecutor` ABC + `FakeRemoteExecutor`

**Files:**
- Create: `mecfs_bio/build_system/wf/remote_executor/base_remote_executor.py`
- Create: `mecfs_bio/build_system/wf/remote_executor/fake_remote_executor.py`
- Test: `test_mecfs_bio/unit/build_system/wf/remote_executor/test_fake_remote_executor.py`

**Interfaces:**
- Consumes: `RemoteJob` (Task 2).
- Produces: `RemoteExecutor.run(self, job: RemoteJob, local_output_dir: Path) -> None` (ABC); `FakeRemoteExecutor(stub_outputs: Mapping[PurePath, str] | None = None)` with `.last_job: RemoteJob | None`, which on `run` records the job and writes each declared `output_file` into `local_output_dir` with stub content (so a consuming Task's asset post-init passes). This behavior is real logic that consumer tests depend on, so it is worth one test.

- [ ] **Step 1: Write failing test.**

```python
from pathlib import Path, PurePath

from mecfs_bio.build_system.wf.remote_executor.fake_remote_executor import FakeRemoteExecutor
from mecfs_bio.build_system.wf.remote_executor.remote_job import RemoteJob, RemoteResources


def test_fake_records_job_and_writes_declared_outputs(tmp_path: Path) -> None:
    ex = FakeRemoteExecutor(stub_outputs={PurePath("out/pip.txt"): "SNP\tPIP\n"})
    job = RemoteJob(
        image="i", commands=["gctb"], input_files={}, s3_inputs={},
        output_files=[PurePath("out/pip.txt")], resources=RemoteResources(1, 1, 1),
    )
    ex.run(job, tmp_path)
    assert ex.last_job is job
    assert (tmp_path / "out/pip.txt").read_text().startswith("SNP")
```

- [ ] **Step 2: Run → FAIL.** `pixi r pytest test_mecfs_bio/unit/build_system/wf/remote_executor/test_fake_remote_executor.py -v`
- [ ] **Step 3: Implement.**

```python
# base_remote_executor.py
from abc import ABC, abstractmethod
from pathlib import Path

from mecfs_bio.build_system.wf.remote_executor.remote_job import RemoteJob


class RemoteExecutor(ABC):
    """Runs a RemoteJob's container commands somewhere and retrieves its outputs."""

    @abstractmethod
    def run(self, job: RemoteJob, local_output_dir: Path) -> None: ...
```

```python
# fake_remote_executor.py
from collections.abc import Mapping
from pathlib import Path, PurePath

from mecfs_bio.build_system.wf.remote_executor.base_remote_executor import RemoteExecutor
from mecfs_bio.build_system.wf.remote_executor.remote_job import RemoteJob


class FakeRemoteExecutor(RemoteExecutor):
    def __init__(self, stub_outputs: Mapping[PurePath, str] | None = None) -> None:
        self._stub_outputs: dict[PurePath, str] = dict(stub_outputs or {})
        self.last_job: RemoteJob | None = None

    def run(self, job: RemoteJob, local_output_dir: Path) -> None:
        self.last_job = job
        for out in job.output_files:
            dest = local_output_dir / out
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_text(self._stub_outputs.get(out, ""))
```

- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Run green** → `... | tee experiments/claude/logs/green_task3.log`; confirm pass. Commit. `git commit -m "feat(remote): RemoteExecutor ABC + FakeRemoteExecutor"`

---

## Task 4: `LocalDockerRemoteExecutor` (system-test seam)

**Files:**
- Create: `mecfs_bio/build_system/wf/remote_executor/local_docker_remote_executor.py`
- Test: `test_mecfs_bio/unit/build_system/wf/remote_executor/test_local_docker_remote_executor.py`

**Interfaces:**
- Consumes: `RemoteJob`, `RemoteExecutor`, `execute_command`.
- Produces: `LocalDockerRemoteExecutor(runner: Callable[[list[str]], str] = execute_command)`; runs the container locally (mounts a temp workdir at `/work`, copies `input_files` in and `output_files` back to `local_output_dir`; `s3_inputs` fetched via AWS CLI, skipped when empty; ignores `resources`).

- [ ] **Step 1: Write failing test** — inject a call-recording runner and assert a `docker run` command was issued. (Real output copy-back is covered by the Task 12 system test, so the unit test stays simple.)

```python
from pathlib import Path, PurePath

from mecfs_bio.build_system.wf.remote_executor.local_docker_remote_executor import LocalDockerRemoteExecutor
from mecfs_bio.build_system.wf.remote_executor.remote_job import RemoteJob, RemoteResources


def test_issues_a_docker_run_command(tmp_path: Path) -> None:
    calls: list[list[str]] = []

    def fake_runner(cmd: list[str]) -> str:
        calls.append(cmd)
        return ""

    ex = LocalDockerRemoteExecutor(runner=fake_runner)
    job = RemoteJob(
        image="busybox:latest", commands=["echo hi"], input_files={}, s3_inputs={},
        output_files=[PurePath("out/pip.txt")], resources=RemoteResources(1, 1, 1),
    )
    # output copy-back would fail because the fake runner produces nothing; that path is
    # exercised for real in the system test, so allow the missing-output error here.
    try:
        ex.run(job, tmp_path)
    except AssertionError:
        pass
    assert any("docker run" in " ".join(c) for c in calls)
```

- [ ] **Step 2: Run → FAIL.**
- [ ] **Step 3: Implement.** Use a `tempfile.TemporaryDirectory` host workdir; `shutil.copy` each `input_files` local→`hostwork/<remote>`; build `docker run --rm -v <hostwork>:/work -w /work <image> bash -lc "<cmd1 && cmd2 ...>"` and pass to `self._runner`; then copy each `output_file` from `hostwork` to `local_output_dir`, asserting existence with a clear message. All params/locals annotated.
- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Run green** → `... | tee experiments/claude/logs/green_task4.log`; confirm pass. Commit. `git commit -m "feat(remote): LocalDockerRemoteExecutor"`

---

## Task 5: `SkyPilotRemoteExecutor` (production AWS)

**Files:**
- Create: `mecfs_bio/build_system/wf/remote_executor/skypilot_remote_executor.py`
- Modify: `pyproject.toml` (add a pinned `skypilot[aws]` to the analysis env; verify it resolves under pixi)
- Test: `test_mecfs_bio/unit/build_system/wf/remote_executor/test_skypilot_remote_executor.py`

**Interfaces:**
- Consumes: `RemoteJob`, `RemoteExecutor`, the SkyPilot SDK (`sky.Task`, `sky.Resources`, `sky.launch`, `sky.get`, `sky.tail_logs`, `sky.down`).
- Produces: `SkyPilotRemoteExecutor(confirm: Callable[[str], bool] = _prompt_confirm, idle_minutes_to_autostop: int = 15)`; pure helpers `build_sky_task(job: RemoteJob) -> "sky.Task"` and `estimate_cost_usd(job: RemoteJob, hours: float) -> float`, unit-tested without AWS. `import sky` goes at module top (SkyPilot is a declared dep); the module is imported lazily by the WF default (Task 7) so unrelated build runs don't pay for it.

- [ ] **Step 1: Write failing test for the pure helper** (no cloud). Define the resource numbers once as locals and assert the built task echoes them:

```python
from pathlib import Path, PurePath

from mecfs_bio.build_system.wf.remote_executor.remote_job import RemoteJob, RemoteResources
from mecfs_bio.build_system.wf.remote_executor.skypilot_remote_executor import build_sky_task


def test_build_sky_task_reflects_resources_and_commands(tmp_path: Path) -> None:
    memory_gb = 192
    vcpus = 24
    disk_gb = 500
    region = "us-east-1"
    job = RemoteJob(
        image="img:1", commands=["gctb --gwfm RC ..."],
        input_files={tmp_path / "x.ma": PurePath("work/x.ma")},
        s3_inputs={"s3://b/Imputed13M/v1": PurePath("work/ref")},
        output_files=[PurePath("work/out")],
        resources=RemoteResources(memory_gb=memory_gb, vcpus=vcpus, disk_gb=disk_gb, region=region),
    )
    task = build_sky_task(job)
    res = list(task.resources)[0]
    assert res.memory == memory_gb
    assert res.cpus == vcpus
    assert res.disk_size == disk_gb
    assert res.region == region
    assert "gctb --gwfm" in task.run
    assert "aws s3" in task.setup
```

- [ ] **Step 2: Run → FAIL.**
- [ ] **Step 3: Implement** (all imports at top, including `import sky`). `build_sky_task` maps `input_files` → `sky.Task.file_mounts`, builds `setup` (an `aws s3 cp --recursive` per `s3_inputs` entry + `docker pull` of `job.image`) and `run` (`docker run --rm -v $(pwd):/work -w /work <image> bash -lc "<commands joined by &&>"`), and sets `sky.Resources(cloud=sky.AWS(), cpus=..., memory=..., disk_size=..., region=...)`. `run()`:

```python
def run(self, job: RemoteJob, local_output_dir: Path) -> None:
    hours_est = 14.0
    prompt = (
        f"Launch {job.resources.vcpus}vCPU/{job.resources.memory_gb}GB on-demand "
        f"(~${estimate_cost_usd(job, hours_est):.0f})? [y/N] "
    )
    if not self._confirm(prompt):
        raise RuntimeError("Remote launch declined by user")
    cluster = f"gwfm-{uuid4().hex[:8]}"
    task = build_sky_task(job)
    try:
        request_id = sky.launch(
            task, cluster_name=cluster,
            idle_minutes_to_autostop=self._idle_minutes_to_autostop, down=True,
        )
        job_id, _ = sky.get(request_id)
        sky.tail_logs(cluster, job_id, follow=True)  # stream ~13h locally
        self._retrieve_outputs(cluster, job.output_files, local_output_dir)
    finally:
        sky.down(cluster)  # guaranteed teardown even on exception/Ctrl-C
```

`_prompt_confirm(prompt: str) -> bool` returns True on a "y" stdin answer OR when `os.environ.get("GWFM_ASSUME_YES") == "1"` (non-interactive override). Module docstring documents the AWS-credential-chain assumption and `GWFM_ASSUME_YES`. Confirm the installed SkyPilot's output-retrieval call for `_retrieve_outputs` (e.g. `sky.Storage`/rsync/scp); if none is direct, have the `run` phase `aws s3 cp` outputs to a run-scoped prefix and download from there.

- [ ] **Step 4: Run helper test → PASS.** The live launch path is validated by a one-time manual maintainer smoke test (Task 12 covers the local path); it is not in CI.
- [ ] **Step 5: Run green** → `... | tee experiments/claude/logs/green_task5.log`; confirm pass. Commit. `git commit -m "feat(remote): SkyPilotRemoteExecutor + pinned skypilot dep"`

---

## Task 6: `ObjectStore` capability (S3 + fake)

**Files:**
- Create: `mecfs_bio/build_system/wf/object_store/base_object_store.py`
- Create: `mecfs_bio/build_system/wf/object_store/s3_object_store.py`
- Create: `mecfs_bio/build_system/wf/object_store/fake_object_store.py`
- Modify: `pyproject.toml` (ensure `boto3` present in analysis env)
- Test: `test_mecfs_bio/unit/build_system/wf/object_store/test_fake_object_store.py`

**Interfaces:**
- Produces: `ObjectHead(size_bytes: int, sha256: str | None)` frozen; `ObjectStore` ABC with `head(self, uri: str) -> ObjectHead | None` and `upload_from_url(self, source_url: str, uri: str) -> str`; `S3ObjectStore(client=...)`; `FakeObjectStore(objects: dict[str, ObjectHead] | None = None)` with `.uploaded: list[tuple[str, str]]`.

- [ ] **Step 1: Write failing test for the fake** (its dedup logic is what Task 8 relies on, so it is worth testing):

```python
from mecfs_bio.build_system.wf.object_store.base_object_store import ObjectHead
from mecfs_bio.build_system.wf.object_store.fake_object_store import FakeObjectStore


def test_fake_head_and_upload_roundtrip() -> None:
    store = FakeObjectStore(objects={"s3://b/present": ObjectHead(size_bytes=10, sha256="a" * 64)})
    assert store.head("s3://b/absent") is None
    assert store.head("s3://b/present").size_bytes == 10
    returned = store.upload_from_url("https://x/file", "s3://b/absent")
    assert ("https://x/file", "s3://b/absent") in store.uploaded
    assert store.head("s3://b/absent") is not None
    assert isinstance(returned, str)
```

- [ ] **Step 2: Run → FAIL.**
- [ ] **Step 3: Implement** ABC + `ObjectHead`, `FakeObjectStore` (in-memory; `upload_from_url` records the call, registers a synthetic `ObjectHead`, returns a sha256 string), and `S3ObjectStore`: `head` via `head_object` (+ `get_object_attributes(ObjectAttributes=["Checksum"])` for the stored SHA-256, `None` if absent); `upload_from_url` streams the URL body via `upload_fileobj(..., ExtraArgs={"StorageClass": "ONEZONE_IA", "ChecksumAlgorithm": "SHA256"})` and returns the stored checksum. `S3ObjectStore` has no unit test (only reachable through the fake + the manual first-staging run) — document that.
- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Run green** → `... | tee experiments/claude/logs/green_task6.log`; confirm pass. Commit. `git commit -m "feat(objectstore): ObjectStore capability with S3 + fake"`

---

## Task 7: Wire capabilities into `WF`, `make_wf`, runner, config

No dedicated unit test (asserting an injected object is returned is cargo-cult); verified by `invoke green` (import-linter + typecheck) and the downstream tests that call `make_wf(remote_executor=fake)`.

**Files:**
- Modify: `mecfs_bio/build_system/wf/base_wf.py`
- Modify: `mecfs_bio/build_system/runner/simple_runner.py:115`
- Modify: `mecfs_bio/analysis/runner/default_runner.py`

**Interfaces:**
- Produces: `WF` gains `remote_executor: RemoteExecutor` and `object_store: ObjectStore`; `make_wf(..., remote_executor: RemoteExecutor | None = None, object_store: ObjectStore | None = None) -> WF` defaults them; new config keys `remote_region`, `remote_s3_bucket`, `gctb_image` read by `default_runner`.

- [ ] **Step 1: Implement.** Add the two attributes to the frozen `WF` and the two params to `make_wf`. To avoid importing SkyPilot on every build, default `remote_executor` via a small lazy factory (import `skypilot_remote_executor` inside the factory) rather than constructing it at import time; default `object_store` to `S3ObjectStore()`. In `simple_runner.py` pass both through `make_wf(...)`. In `default_runner.py` add `_REMOTE_REGION_KEY` / `_REMOTE_BUCKET_KEY` / `_GCTB_IMAGE_KEY` readers mirroring `_get_asset_root_path`, and update its module docstring listing the new optional config keys.
- [ ] **Step 2: Run green** → `... | tee experiments/claude/logs/green_task7.log`; confirm pass (watch the import-linter — keep the new `wf` submodules within their layer).
- [ ] **Step 3: Commit.** `git commit -m "feat(wf): wire remote_executor + object_store capabilities"`

---

## Task 8: `StageGwfmReferenceTask`

**Files:**
- Create: `mecfs_bio/build_system/task/sbayesrc/stage_gwfm_reference_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/sbayesrc/test_stage_gwfm_reference_task.py`

**Interfaces:**
- Consumes: `GWFM_REFERENCE_BUNDLE`, `GWFM_REFERENCE_VERSION`, the `MARKER_*_KEY` constants (Task 1); `ObjectStore`/`ObjectHead` (Task 6); `FileAsset`, `Fetch`, `WF`.
- Produces: `StageGwfmReferenceTask(bucket: str)`; `execute` returns a `FileAsset` marker whose JSON keys are exactly `MARKER_VERSION_KEY`/`MARKER_PREFIX_KEY`/`MARKER_FILES_KEY` (deterministic across machines from the pinned bundle).

- [ ] **Step 1: Write failing tests** (inject a `FakeObjectStore` through `wf`). Import `MARKER_VERSION_KEY` and use the SAME constant in the assertion, so the reader sees where the key comes from:

```python
from pathlib import Path

from mecfs_bio.build_system.task.sbayesrc.gctb_gwfm_constants import MARKER_VERSION_KEY
# ... construct wf whose object_store is a FakeObjectStore pre-populated with all bundle
# files at their expected size ...


def test_skips_upload_when_bucket_already_populated(tmp_path: Path, wf_all_present) -> None:
    task = StageGwfmReferenceTask(bucket="mybucket")
    asset = task.execute(tmp_path, fetch=_noop_fetch, wf=wf_all_present)
    assert wf_all_present.object_store.uploaded == []
    assert MARKER_VERSION_KEY in Path(asset.path).read_text()


def test_uploads_only_missing_files(tmp_path: Path, wf_one_missing) -> None:
    task = StageGwfmReferenceTask(bucket="mybucket")
    task.execute(tmp_path, fetch=_noop_fetch, wf=wf_one_missing)
    assert len(wf_one_missing.object_store.uploaded) == 1
```

- [ ] **Step 2: Run → FAIL.**
- [ ] **Step 3: Implement.** For each `ReferenceBundleFile`: `uri = f"s3://{bucket}/sbayesrc/ld/{GWFM_REFERENCE_VERSION}/{f.filename}"`; `head = wf.object_store.head(uri)`; upload (via `upload_from_url`) iff `head is None` or `head.size_bytes != f.size_bytes` or (`f.sha256` and `head.sha256 != f.sha256`); log each returned sha256 so it can be pasted back into the constant. Build the marker dict with the `MARKER_*_KEY` constants (`json.dumps(..., sort_keys=True)`), write it, return `FileAsset`. Assert the bundle is non-empty (fail fast).
- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Run green** → `... | tee experiments/claude/logs/green_task8.log`; confirm pass. Commit. `git commit -m "feat(gwfm): StageGwfmReferenceTask with S3 dedup + deterministic marker"`

> After the FIRST real staging run, paste the logged per-file SHA-256s into `GWFM_REFERENCE_BUNDLE` and commit, so future traces are checksum-verified.

---

## Task 9: `SumstatsToCojoMaTask`

**Files:**
- Create: `mecfs_bio/build_system/task/sbayesrc/sumstats_to_cojo_ma_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/sbayesrc/test_sumstats_to_cojo_ma_task.py`

**Interfaces:**
- Consumes: a sumstats parquet dep (read via `scan_dataframe_asset`); `Fetch`, `WF`, `FileAsset`.
- Produces: `SumstatsToCojoMaTask(...)` writing a `.ma` with GCTB/COJO columns `SNP A1 A2 freq b se p N` (tab-separated, header). Trait/project derived from the dep's meta via an isinstance assertion (not passed as args).

- [ ] **Step 1: Write failing test** on a tiny 3-row polars frame → read the `.ma` back and assert the 8 columns in order with correct values.
- [ ] **Step 2: Run → FAIL.**
- [ ] **Step 3: Implement** with polars: fetch dep, select/rename to the 8 COJO columns using column-name constants (add any missing to `constants/gwaslab_constants.py`), `write_csv(separator="\t")`. Return `FileAsset`. All params annotated.
- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Run green** → `... | tee experiments/claude/logs/green_task9.log`; confirm pass. Commit. `git commit -m "feat(gwfm): SumstatsToCojoMaTask (parquet -> .ma)"`

---

## Task 10: `GctbFineMapTask`

**Files:**
- Create: `mecfs_bio/build_system/task/sbayesrc/gctb_fine_map_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/sbayesrc/test_gctb_fine_map_task.py`

**Interfaces:**
- Consumes: `SumstatsToCojoMaTask` (`.ma`), `StageGwfmReferenceTask` (marker); CLI templates + resource constants (Task 1); `WF.remote_executor`; `DirectoryAsset`.
- Produces: `GctbFineMapTask(ma_task: Task, reference_task: Task, threads: int = DEFAULT_VCPUS)`; `deps == [ma_task, reference_task]`; `execute` builds a `RemoteJob` (commands from templates; `input_files` = the `.ma`; `s3_inputs` = the marker's prefix → `PurePath("work/ref")`; `output_files = [PurePath("work/out")]`; `resources = RemoteResources(DEFAULT_MEMORY_GB, threads, DEFAULT_DISK_GB, region)`), calls `wf.remote_executor.run(job, scratch_dir)`, returns `DirectoryAsset(scratch_dir / "work/out")`. Whether `--make-ldm-eigen` is emitted follows the Task 1 finding (authored in, not a runtime flag).

- [ ] **Step 1: Write failing test** with `FakeRemoteExecutor` (stub the `out/` contents so `DirectoryAsset` passes). Assert on the essential invariants only — not the exact command count:

```python
from pathlib import Path, PurePath

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.wf.remote_executor.fake_remote_executor import FakeRemoteExecutor
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.build_system.task.sbayesrc import gctb_gwfm_constants as c


def test_dispatches_a_wellformed_gwfm_job(tmp_path: Path, ma_asset, marker_asset) -> None:
    fake = FakeRemoteExecutor(stub_outputs={PurePath("work/out/snpRes.txt"): "x"})
    wf = make_wf(remote_executor=fake)
    task = GctbFineMapTask(ma_task=_const(ma_asset), reference_task=_const(marker_asset))
    asset = task.execute(tmp_path, fetch=_fetch_for(ma_asset, marker_asset), wf=wf)

    job = fake.last_job
    assert any("--gwfm RC" in cmd for cmd in job.commands)  # essential step present
    assert any(str(remote).endswith(".ma") for remote in job.input_files.values())
    assert job.resources.memory_gb == c.DEFAULT_MEMORY_GB
    assert isinstance(asset, DirectoryAsset)
```

- [ ] **Step 2: Run → FAIL.**
- [ ] **Step 3: Implement** (read the marker JSON via `fetch` using the `MARKER_*_KEY` constants, format templates with concrete `work/`-relative paths, assemble the `RemoteJob`, call the executor, return the `DirectoryAsset`).
- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Run green** → `... | tee experiments/claude/logs/green_task10.log`; confirm pass. Commit. `git commit -m "feat(gwfm): GctbFineMapTask dispatches remote gctb job"`

---

## Task 11: GCTB Docker image + publish task

**Files:**
- Create: `docker/gctb/Dockerfile`
- Create: `docker/gctb/LICENSE-GCTB` (copy of GCTB's MIT license + copyright)
- Modify: `tasks.py` (add `build_push_gctb_image`)
- Test: `test_mecfs_bio/system/test_gctb_image.py`

**Interfaces:**
- Consumes: `GCTB_VERSION`, `GCTB_BINARY_URL`, `GCTB_BINARY_SHA256` (Task 1); `execute_command`.
- Produces: a runnable image tag `<registry>/gctb:<version>`; `invoke build-push-gctb-image` (maintainer-only).

- [ ] **Step 1: Write the Dockerfile.**

```dockerfile
FROM debian:stable-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgomp1 ca-certificates curl unzip \
    && rm -rf /var/lib/apt/lists/*
ARG GCTB_URL
ARG GCTB_SHA256
RUN curl -fsSL "$GCTB_URL" -o /tmp/gctb.zip \
    && echo "$GCTB_SHA256  /tmp/gctb.zip" | sha256sum -c - \
    && unzip /tmp/gctb.zip -d /opt/gctb && rm /tmp/gctb.zip \
    && ln -s "$(find /opt/gctb -name gctb -type f | head -1)" /usr/local/bin/gctb
COPY LICENSE-GCTB /opt/gctb/LICENSE-GCTB
ENTRYPOINT []
```

- [ ] **Step 2: Write the failing system test** using `execute_command` (not `subprocess`): build the image, run `gctb` in it, assert the returned output contains the GCTB banner.

```python
from mecfs_bio.build_system.task.sbayesrc import gctb_gwfm_constants as c
from mecfs_bio.util.subproc.run_command import execute_command


def test_gctb_image_builds_and_runs() -> None:
    execute_command([
        "docker", "build",
        "--build-arg", f"GCTB_URL={c.GCTB_BINARY_URL}",
        "--build-arg", f"GCTB_SHA256={c.GCTB_BINARY_SHA256}",
        "-t", "gctb:test", "docker/gctb",
    ])
    output = execute_command(["docker", "run", "--rm", "gctb:test", "gctb"])
    assert "GCTB" in output
```

- [ ] **Step 3: Run → FAIL, then implement** the `build_push_gctb_image` invoke task (builds with the two build-args via `execute_command`, tags `<registry>/gctb:<GCTB_VERSION>` from config, pushes). Document it as maintainer-only.
- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Run green** → `... | tee experiments/claude/logs/green_task11.log`; confirm pass. Commit. `git commit -m "feat(gwfm): minimal MIT-compliant GCTB image + publish task"`

---

## Task 12: End-to-end system test via `LocalDockerRemoteExecutor`, with a known-truth toy model

Tests not just that the image runs, but that GWFM recovers planted causal SNPs.

**Files:**
- Create: `test_mecfs_bio/system/test_gctb_gwfm.py`
- Create: `test_mecfs_bio/system/test_data/gwfm_toy/` (generator script + generated small inputs)

**Interfaces:**
- Consumes: `GctbFineMapTask`, `LocalDockerRemoteExecutor`, the `gctb:test` image, toy data.

- [ ] **Step 1: Write a toy-data generator** (`experiments/claude/` script, outputs committed under `test_data/gwfm_toy/`): synthesize a small genotype matrix for one LD block (a few hundred SNPs), pick `k` truly causal SNPs with defined effect sizes, compute the LD matrix `R` and marginal GWAS stats consistent with `b_marginal = R @ b_true` plus small noise at a chosen `N`, and write the `.ma`, LD reference, annotation, and gene-map in the formats `gctb` expects (confirm formats against the Task 1 findings). Sized so `gctb --gwfm RC` finishes in seconds.
- [ ] **Step 2: Write the test** — run `GctbFineMapTask.execute` with `make_wf(remote_executor=LocalDockerRemoteExecutor())` (the toy reference is passed as local `input_files`, so `s3_inputs` is empty); parse the PIP output and assert the planted causal SNPs receive high PIP:

```python
CAUSAL_SNPS = ["rs_causal_1", "rs_causal_2"]
PIP_THRESHOLD = 0.5


def test_gwfm_recovers_planted_causal_snps(tmp_path: Path) -> None:
    asset = _run_gwfm_on_toy(tmp_path)  # builds task + local-docker wf, returns DirectoryAsset
    pip = _read_pip_table(asset.path)   # polars: SNP -> PIP
    for snp in CAUSAL_SNPS:
        assert pip[snp] > PIP_THRESHOLD
```

- [ ] **Step 3: Run → PASS** (`pixi r pytest test_mecfs_bio/system/test_gctb_gwfm.py -v`). Tune `N`/effect sizes if PIPs are marginal, keeping runtime in seconds.
- [ ] **Step 4: Run green** → `... | tee experiments/claude/logs/green_task12.log`; confirm pass. Commit. `git commit -m "test(gwfm): end-to-end GWFM recovers planted causal SNPs on toy data"`

---

## Task 13: Collaborator setup notes + config template

**Files:**
- Create: `default_runner_config.example.yaml` (committed; `default_runner_config.yaml` itself stays gitignored)
- Create: `experiments/claude/gwfm_collaborator_setup.md` (NOT under `docs/` — that is the human-written published site)

- [ ] **Step 1** Write the example config with the remote-exec keys commented (`remote_region`, `remote_s3_bucket`, `gctb_image`), and the setup notes: `pixi install`, `aws configure`/SSO, `sky check`, fill the gitignored config, run. Note `GWFM_ASSUME_YES=1` for non-interactive runs; credentials come only from the standard AWS chain, never the repo.
- [ ] **Step 2: Run green** → `... | tee experiments/claude/logs/green_task13.log`; confirm pass. Commit. `git commit -m "docs(gwfm): collaborator setup notes + example remote config (experiments/claude)"`

---

## Self-Review

**Spec coverage:** Approach-1 WF seam → Tasks 2–7,10; SkyPilot/AWS → Task 5; public MIT image → Task 11; reference bundle (LD+eigen+annot+gene-map) staging with deterministic marker + S3 checksum dedup → Task 8; `.ma` conversion → Task 9; on-demand-only / no `use_spot` → Task 2 + Global Constraints; cost/failure safety (teardown finally, idle-autostop, pre-launch confirm, streamed logs) → Task 5; One Zone-IA same-region → Tasks 6/Constraints; testing (fake unit + local-docker system with known-truth PIP check) → Tasks 3,4,10,12; config/credentials (AWS chain, gitignored config) → Tasks 7,13; future local-annotation asset → out of scope (recorded in spec). Recon-dependent externals → Task 1.

**Placeholder scan:** deferred values (bundle SHA-256s, exact gctb flags, disk size, `--make-ldm-eigen` inclusion, SkyPilot output-retrieval call, toy-input formats) are explicitly produced by Task 1 / recorded after first staging, and referenced by name — not left vague in consumer tasks.

**Type consistency:** `RemoteJob`/`RemoteResources` fields, `RemoteExecutor.run(job, local_output_dir)`, `ObjectStore.head/upload_from_url`, `ObjectHead(size_bytes, sha256)`, and the marker JSON keys (`MARKER_VERSION_KEY`/`MARKER_PREFIX_KEY`/`MARKER_FILES_KEY`) are used identically across Tasks 2–12.

## Findings recorded outside the repo (Task 1 → `experiments/claude/`)
1. GCTB license → **MIT, redistribution OK** (already resolved; image ships `LICENSE-GCTB`).
2. GCTB resume support → informational; gates any future spot work.
3. Precomputed `eigen/` usable directly → decides whether Task 10 emits `--make-ldm-eigen`; plus VM disk sizing.
4. Exact GWFM annotation + gene-map files/sizes on gctbhub, and `gctb` input formats for the toy generator.
5. Current S3 One Zone-IA pricing for the chosen region.
6. Exact SkyPilot output-retrieval API for the installed version (Task 5 Step 3).
