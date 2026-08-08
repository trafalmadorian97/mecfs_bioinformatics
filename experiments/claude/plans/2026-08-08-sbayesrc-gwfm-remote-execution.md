# SBayesRC GWFM Remote Execution — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run SBayesRC genome-wide fine-mapping (GWFM) for a trait on a transient AWS instance while the build graph runs locally, with results flowing back as a normal cached build asset.

**Architecture:** A new `RemoteExecutor` capability on the build system's `WF` object dispatches a self-contained container command (the `gctb` GWFM run) to a SkyPilot-provisioned instance; small per-run inputs are shipped from the laptop, the large shared reference bundle is pulled from S3 on the instance, and outputs are retrieved into the local scratch dir so the verifying-trace rebuilder caches them normally. The fine-mapping Task only builds a `RemoteJob` — it is agnostic to the provisioning mechanism.

**Tech Stack:** Python 3.13, attrs (`@frozen`), polars, boto3 (S3), SkyPilot (AWS), Docker, pytest, pixi task runner, structlog.

## Global Constraints

- All commands run via `pixi r <command>`; Python via `pixi r python <script>`. (CLAUDE.md)
- After any significant change run `pixi r invoke green` (lint, format, spellcheck, link-check, import-check, typecheck, test); capture output to a logfile — testmon can report "no tests ran". (project memory)
- In docstrings: no backticks around inline code, no RST. (CLAUDE.md)
- Prefer `Path` for filesystem paths, `PurePath` for base-relative paths; convert to `str` only at serialization boundaries. (project memory)
- Prefer polars over pandas for new dataframe code. (project memory)
- Return `@frozen` objects with named attributes, not bare tuples. (project memory)
- Never monkeypatch/mock — inject dependencies (executor, sleep, object store, subprocess runner) as params with production defaults; tests pass their own. (project memory)
- Enforce cross-field invariants in `__attrs_post_init__`; assert and fail fast on invalid input. (project memory)
- No `skipif` on library presence; no assertions on error-message/log text; test at the Task level where practical. (project memory)
- Type small fixed-value string params as `Literal` aliases. (project memory)
- On-demand instances only. Do NOT add spot support unless Task 1 finds explicit GCTB MCMC-resumption evidence. `use_spot` must not appear in any interface.
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

GWFM domain:
- `mecfs_bio/build_system/task/sbayesrc/gctb_gwfm_constants.py` — pinned versions, URLs, bundle manifest, CLI templates, resource/disk sizing.
- `mecfs_bio/build_system/task/sbayesrc/stage_gwfm_reference_task.py` — `StageGwfmReferenceTask`.
- `mecfs_bio/build_system/task/sbayesrc/sumstats_to_cojo_ma_task.py` — `SumstatsToCojoMaTask`.
- `mecfs_bio/build_system/task/sbayesrc/gctb_fine_map_task.py` — `GctbFineMapTask`.
- `mecfs_bio/build_system/meta/gwfm_reference_marker_meta.py` — marker file meta (or reuse `ReferenceFileMeta`; Task 1 decides).

Wiring / config / infra:
- `mecfs_bio/build_system/wf/base_wf.py` — add `remote_executor`, `object_store` to `WF` + `make_wf`.
- `mecfs_bio/build_system/runner/simple_runner.py` — pass the two capabilities into `make_wf`.
- `mecfs_bio/analysis/runner/default_runner.py` — read remote-exec config (region, bucket, image, non-interactive flag).
- `docker/gctb/Dockerfile` — minimal GCTB image.
- `tasks.py` — `build_push_gctb_image` invoke task.

Tests mirror under `test_mecfs_bio/unit/build_system/...` and one `test_mecfs_bio/system/test_gctb_gwfm.py`.

---

## Task 1: Reconnaissance — pin external constants

**Files:**
- Create: `mecfs_bio/build_system/task/sbayesrc/gctb_gwfm_constants.py`
- Create: `experiments/claude/design_specs/2026-08-08-gwfm-recon-findings.md`
- Test: `test_mecfs_bio/unit/build_system/task/sbayesrc/test_gctb_gwfm_constants.py`

**Interfaces:**
- Produces: module-level constants consumed by every later GWFM task —
  `GCTB_VERSION: str`, `GCTB_BINARY_URL: str`, `GCTB_BINARY_SHA256: str`,
  `GWFM_REFERENCE_VERSION: str` (e.g. `"Imputed13M/v1"`),
  `ReferenceBundleFile` (frozen: `filename: str`, `source_url: str`, `size_bytes: int`, `sha256: str | None`),
  `GWFM_REFERENCE_BUNDLE: tuple[ReferenceBundleFile, ...]`,
  `GCTB_MAKE_LDM_EIGEN_TEMPLATE`, `GCTB_GWFM_TEMPLATE`, `GCTB_CS_TEMPLATE` (str templates with named `{}` fields),
  `USES_PRECOMPUTED_EIGEN: bool`, `DEFAULT_MEMORY_GB: int`, `DEFAULT_VCPUS: int`, `DEFAULT_DISK_GB: int`,
  `GCTB_SUPPORTS_RESUME: bool`.

- [ ] **Step 1: Investigate and record findings.** Confirm, and write into the recon-findings doc with a source link for each:
  1. Exact gctbhub GWFM reference files + directory-listing URLs under `https://gctbhub.cloud.edu.au/data/SBayesRC/resources/GWFM/LD/Imputed13M/` and the annotation + gene-map location (check `.../resources/v2.0/Annotation/` and the GWFM tutorial). Record each file's `Content-Length` (HEAD) as `size_bytes`.
  2. The precompiled binary URL `https://gctbhub.cloud.edu.au/software/gctb/download/gctb_2.5.5_Linux.zip` and its sha256 (download once locally, `sha256sum`, record).
  3. Exact GWFM CLI from the tutorial (the three `gctb` calls) and whether `--gwfm RC` can consume the precomputed `eigen/` directly (set `USES_PRECOMPUTED_EIGEN` and adjust whether step-1 `--make-ldm-eigen` is needed).
  4. Whether `gctb` GWFM supports checkpoint/resume (search the GCTB repo/docs for "resume"/".rds"/checkpoint). Set `GCTB_SUPPORTS_RESUME`; if False, spot stays out of scope.
  5. VM disk sizing: reference zip + unzip + eigen intermediates → set `DEFAULT_DISK_GB` (plan ~500).

- [ ] **Step 2: Write the constants module** using the findings. Reference-bundle `sha256` fields start as `None` (filled and committed after the first real staging run in Task 8). Example shape:

```python
from attrs import frozen

GCTB_VERSION = "2.5.5"
GCTB_BINARY_URL = "https://gctbhub.cloud.edu.au/software/gctb/download/gctb_2.5.5_Linux.zip"
GCTB_BINARY_SHA256 = "<recorded>"  # from Step 1.2

GWFM_REFERENCE_VERSION = "Imputed13M/v1"

@frozen
class ReferenceBundleFile:
    filename: str
    source_url: str
    size_bytes: int
    sha256: str | None  # None until first staging records it; then committed

GWFM_REFERENCE_BUNDLE: tuple[ReferenceBundleFile, ...] = (
    ReferenceBundleFile("ukbEUR_13M_FullLDM.zip", "https://.../ukbEUR_13M_FullLDM.zip", 0, None),
    # eigen/ files, annot file, gene-map file ...
)

# {ldm}, {ma}, {annot}, {gene_map}, {out}, {threads}, {pip}, {pep} filled at run time.
GCTB_MAKE_LDM_EIGEN_TEMPLATE = "gctb --ldm {ldm} --gwas-summary {ma} --make-ldm-eigen --thread {threads} --out {out}"
GCTB_GWFM_TEMPLATE = "gctb --gwfm RC --ldm-eigen {ldm} --gwas-summary {ma} --annot {annot} --gene-map {gene_map} --thread {threads} --out {out}"
GCTB_CS_TEMPLATE = "gctb --cs --pwld-file {pwld} --pip {pip} --pep {pep} --gene-map {gene_map} --mcmc-samples {out} --out {out}"

USES_PRECOMPUTED_EIGEN = False   # per Step 1.3
DEFAULT_MEMORY_GB = 192
DEFAULT_VCPUS = 24
DEFAULT_DISK_GB = 500
GCTB_SUPPORTS_RESUME = False     # per Step 1.4
```

- [ ] **Step 3: Write invariant tests.**

```python
import re
from mecfs_bio.build_system.task.sbayesrc import gctb_gwfm_constants as c

def test_binary_checksum_is_hex_sha256():
    assert re.fullmatch(r"[0-9a-f]{64}", c.GCTB_BINARY_SHA256)

def test_bundle_is_nonempty_and_urls_are_https():
    assert c.GWFM_REFERENCE_BUNDLE
    for f in c.GWFM_REFERENCE_BUNDLE:
        assert f.source_url.startswith("https://")
        assert f.size_bytes > 0

def test_cli_templates_have_expected_fields():
    assert "{ma}" in c.GCTB_GWFM_TEMPLATE and "{annot}" in c.GCTB_GWFM_TEMPLATE
```

- [ ] **Step 4: Run tests.** `pixi r pytest test_mecfs_bio/unit/build_system/task/sbayesrc/test_gctb_gwfm_constants.py -v` → PASS.
- [ ] **Step 5: Commit.** `git add mecfs_bio/build_system/task/sbayesrc/gctb_gwfm_constants.py experiments/claude/design_specs/2026-08-08-gwfm-recon-findings.md test_mecfs_bio/unit/build_system/task/sbayesrc/test_gctb_gwfm_constants.py && git commit -m "feat(gwfm): pin GCTB GWFM constants from recon"`

---

## Task 2: `RemoteResources` and `RemoteJob` types

**Files:**
- Create: `mecfs_bio/build_system/wf/remote_executor/remote_job.py`
- Test: `test_mecfs_bio/unit/build_system/wf/remote_executor/test_remote_job.py`

**Interfaces:**
- Produces: `RemoteResources(memory_gb:int, vcpus:int, disk_gb:int, region:str|None=None)`;
  `RemoteJob(image:str, commands:Sequence[str], input_files:Mapping[Path,PurePath], s3_inputs:Mapping[str,PurePath], output_files:Sequence[PurePath], resources:RemoteResources)`.

- [ ] **Step 1: Write failing tests.**

```python
from pathlib import Path, PurePath
import pytest
from mecfs_bio.build_system.wf.remote_executor.remote_job import RemoteResources, RemoteJob

def _res(): return RemoteResources(memory_gb=192, vcpus=24, disk_gb=500)

def test_remote_job_holds_fields():
    job = RemoteJob(image="img@sha256:abc", commands=["gctb --help"],
                    input_files={Path("/tmp/x.ma"): PurePath("work/x.ma")},
                    s3_inputs={"s3://b/k": PurePath("work/ref")},
                    output_files=[PurePath("work/out")], resources=_res())
    assert job.commands[0].startswith("gctb")

def test_positive_resources_enforced():
    with pytest.raises(AssertionError):
        RemoteResources(memory_gb=0, vcpus=24, disk_gb=500)

def test_empty_commands_rejected():
    with pytest.raises(AssertionError):
        RemoteJob(image="i", commands=[], input_files={}, s3_inputs={},
                  output_files=[PurePath("o")], resources=_res())
```

- [ ] **Step 2: Run → FAIL** (module missing). `pixi r pytest test_mecfs_bio/unit/build_system/wf/remote_executor/test_remote_job.py -v`
- [ ] **Step 3: Implement.**

```python
from collections.abc import Mapping, Sequence
from pathlib import Path, PurePath
from attrs import frozen

# NOTE: no `use_spot` / `instance_type` — see design spec. If a named-instance
# override is added later, refactor RemoteJob.resources to a
# `RemoteResources | ExplicitInstance` union rather than adding optional fields.

@frozen
class RemoteResources:
    memory_gb: int
    vcpus: int
    disk_gb: int
    region: str | None = None

    def __attrs_post_init__(self):
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

    def __attrs_post_init__(self):
        assert self.commands, "RemoteJob.commands must be non-empty"
        assert self.output_files, "RemoteJob.output_files must be non-empty"
```

- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Commit.** `git commit -m "feat(remote): RemoteJob/RemoteResources types"`

---

## Task 3: `RemoteExecutor` ABC + `FakeRemoteExecutor`

**Files:**
- Create: `mecfs_bio/build_system/wf/remote_executor/base_remote_executor.py`
- Create: `mecfs_bio/build_system/wf/remote_executor/fake_remote_executor.py`
- Test: `test_mecfs_bio/unit/build_system/wf/remote_executor/test_fake_remote_executor.py`

**Interfaces:**
- Consumes: `RemoteJob` (Task 2).
- Produces: `RemoteExecutor.run(self, job: RemoteJob, local_output_dir: Path) -> None` (ABC);
  `FakeRemoteExecutor(stub_outputs: Mapping[PurePath, str] = {})` with attribute `.last_job: RemoteJob | None`, which on `run` records the job and writes each declared `output_file` into `local_output_dir` (creating parent dirs) with stub content, so a consuming Task's `DirectoryAsset`/`FileAsset` post-init passes.

- [ ] **Step 1: Write failing test.**

```python
from pathlib import Path, PurePath
from mecfs_bio.build_system.wf.remote_executor.fake_remote_executor import FakeRemoteExecutor
from mecfs_bio.build_system.wf.remote_executor.remote_job import RemoteJob, RemoteResources

def test_fake_records_job_and_writes_outputs(tmp_path):
    ex = FakeRemoteExecutor(stub_outputs={PurePath("out/pip.txt"): "SNP\tPIP\n"})
    job = RemoteJob(image="i", commands=["gctb"], input_files={}, s3_inputs={},
                    output_files=[PurePath("out/pip.txt")],
                    resources=RemoteResources(1, 1, 1))
    ex.run(job, tmp_path)
    assert ex.last_job is job
    assert (tmp_path / "out/pip.txt").read_text().startswith("SNP")
```

- [ ] **Step 2: Run → FAIL.**
- [ ] **Step 3: Implement ABC + fake.**

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
    def __init__(self, stub_outputs: Mapping[PurePath, str] | None = None):
        self._stub_outputs = dict(stub_outputs or {})
        self.last_job: RemoteJob | None = None

    def run(self, job: RemoteJob, local_output_dir: Path) -> None:
        self.last_job = job
        for out in job.output_files:
            dest = local_output_dir / out
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_text(self._stub_outputs.get(out, ""))
```

- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Commit.** `git commit -m "feat(remote): RemoteExecutor ABC + FakeRemoteExecutor"`

---

## Task 4: `LocalDockerRemoteExecutor` (system-test seam)

**Files:**
- Create: `mecfs_bio/build_system/wf/remote_executor/local_docker_remote_executor.py`
- Test: `test_mecfs_bio/unit/build_system/wf/remote_executor/test_local_docker_remote_executor.py`

**Interfaces:**
- Consumes: `RemoteJob`, `RemoteExecutor`, `execute_command` (`mecfs_bio.util.subproc.run_command`).
- Produces: `LocalDockerRemoteExecutor(runner: Callable[[list[str]], str] = execute_command)`; runs the container locally, mounting a temp work dir at `/work`, copying `input_files` in and `output_files` back to `local_output_dir`. `s3_inputs` are fetched with the AWS CLI inside/outside the container (skipped when empty). Ignores `resources` (local box).

- [ ] **Step 1: Write failing test** (inject a fake runner; assert the emitted docker command mounts a workdir, references the image, and runs the job commands; assert declared outputs pre-created by the fake runner are copied back).

```python
from pathlib import Path, PurePath
from mecfs_bio.build_system.wf.remote_executor.local_docker_remote_executor import LocalDockerRemoteExecutor
from mecfs_bio.build_system.wf.remote_executor.remote_job import RemoteJob, RemoteResources

def test_builds_docker_command_and_copies_outputs(tmp_path):
    calls = []
    def fake_runner(cmd):
        calls.append(cmd)
        # simulate the container producing its declared output in the mounted workdir
        joined = " ".join(cmd)
        if "docker run" in joined:
            # workdir mount is "<hostwork>:/work"; find it
            mount = [t for t in joined.split() if t.endswith(":/work")][0]
            host = Path(mount.split(":/work")[0])
            (host / "out").mkdir(parents=True, exist_ok=True)
            (host / "out/pip.txt").write_text("ok")
        return ""
    ex = LocalDockerRemoteExecutor(runner=fake_runner)
    job = RemoteJob(image="busybox:latest", commands=["echo hi"], input_files={},
                    s3_inputs={}, output_files=[PurePath("out/pip.txt")],
                    resources=RemoteResources(1, 1, 1))
    ex.run(job, tmp_path)
    assert any("docker run" in " ".join(c) for c in calls)
    assert (tmp_path / "out/pip.txt").read_text() == "ok"
```

- [ ] **Step 2: Run → FAIL.**
- [ ] **Step 3: Implement** (use a `tempfile.TemporaryDirectory` as host workdir; `shutil.copy` inputs to `hostwork/<remote>`; build `docker run --rm -v hostwork:/work -w /work <image> bash -lc "<cmd1 && cmd2 ...>"` via the injected `runner`; then copy each `output_file` from `hostwork` to `local_output_dir`, asserting it exists with a clear message).
- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Commit.** `git commit -m "feat(remote): LocalDockerRemoteExecutor"`

---

## Task 5: `SkyPilotRemoteExecutor` (production AWS)

**Files:**
- Create: `mecfs_bio/build_system/wf/remote_executor/skypilot_remote_executor.py`
- Modify: `pyproject.toml` (add `skypilot[aws]` to the analysis env deps — pin a version; verify it resolves in pixi)
- Test: `test_mecfs_bio/unit/build_system/wf/remote_executor/test_skypilot_remote_executor.py`

**Interfaces:**
- Consumes: `RemoteJob`, `RemoteExecutor`, constants (disk/mem), the SkyPilot SDK (`sky.Task`, `sky.Resources`, `sky.launch`, `sky.get`, `sky.tail_logs`, `sky.down`).
- Produces: `SkyPilotRemoteExecutor(confirm: Callable[[str], bool] = _prompt_confirm, idle_minutes_to_autostop: int = 15)`; plus pure helper `build_sky_task(job: RemoteJob) -> "sky.Task"` and `estimate_cost_usd(job: RemoteJob, hours: float) -> float` that are unit-tested without touching AWS.

- [ ] **Step 1: Write failing tests for the pure helpers only** (no cloud).

```python
from pathlib import Path, PurePath
from mecfs_bio.build_system.wf.remote_executor.remote_job import RemoteJob, RemoteResources
from mecfs_bio.build_system.wf.remote_executor.skypilot_remote_executor import build_sky_task

def test_build_sky_task_sets_resources_and_run(tmp_path):
    job = RemoteJob(image="img:1", commands=["gctb --gwfm RC ..."],
                    input_files={tmp_path / "x.ma": PurePath("work/x.ma")},
                    s3_inputs={"s3://b/Imputed13M/v1": PurePath("work/ref")},
                    output_files=[PurePath("work/out")],
                    resources=RemoteResources(memory_gb=192, vcpus=24, disk_gb=500, region="us-east-1"))
    task = build_sky_task(job)
    res = list(task.resources)[0]
    assert res.memory == 192 and res.cpus == 24 and res.disk_size == 500
    assert res.region == "us-east-1"
    assert "gctb --gwfm" in task.run
    assert "aws s3" in task.setup  # pulls s3_inputs
```

- [ ] **Step 2: Run → FAIL.**
- [ ] **Step 3: Implement.** `build_sky_task` maps `input_files` → `sky.Task.file_mounts` (local→remote), builds a `setup` that `aws s3 cp --recursive` each `s3_inputs` entry and `docker pull`s `job.image`, and a `run` that `docker run --rm -v $(pwd):/work -w /work <image> bash -lc "<commands joined by &&>"`, then declares `resources = sky.Resources(cloud=sky.AWS(), cpus=job.resources.vcpus, memory=job.resources.memory_gb, disk_size=job.resources.disk_gb, region=job.resources.region, image_id=None)`. `run()`:

```python
def run(self, job, local_output_dir):
    hours_est = 14.0
    if not self.confirm(f"Launch {job.resources.vcpus}vCPU/{job.resources.memory_gb}GB "
                        f"on-demand (~${estimate_cost_usd(job, hours_est):.0f})? [y/N] "):
        raise RuntimeError("Remote launch declined by user")
    import sky
    cluster = f"gwfm-{uuid4().hex[:8]}"
    task = build_sky_task(job)
    try:
        req = sky.launch(task, cluster_name=cluster,
                         idle_minutes_to_autostop=self.idle_minutes_to_autostop, down=True)
        job_id, _ = sky.get(req)
        sky.tail_logs(cluster, job_id, follow=True)   # stream ~13h locally
        _rsync_outputs_back(cluster, job.output_files, local_output_dir)  # sky rsync/scp helper
    finally:
        sky.down(cluster)   # guaranteed teardown even on exception/Ctrl-C
```

`_prompt_confirm` returns True when stdin says yes OR when env var `GWFM_ASSUME_YES=1` (non-interactive override). Add a module docstring documenting the AWS-credential-chain assumption and the `GWFM_ASSUME_YES` escape hatch. Verify the exact SkyPilot output-retrieval call (`sky.Storage`/`rsync_down`/`scp`) against the installed version and implement `_rsync_outputs_back` accordingly; if SkyPilot lacks a direct download, have the `run` phase `aws s3 cp` outputs to a run-scoped S3 prefix and download from there.

- [ ] **Step 4: Run helper tests → PASS.** Live launch is validated manually (Task 12 covers the local path; a real cloud smoke test is run once by the maintainer, not in CI).
- [ ] **Step 5: Commit.** `git commit -m "feat(remote): SkyPilotRemoteExecutor + pinned skypilot dep"`

---

## Task 6: `ObjectStore` capability (S3 + fake)

**Files:**
- Create: `mecfs_bio/build_system/wf/object_store/base_object_store.py`
- Create: `mecfs_bio/build_system/wf/object_store/s3_object_store.py`
- Create: `mecfs_bio/build_system/wf/object_store/fake_object_store.py`
- Modify: `pyproject.toml` (ensure `boto3` present in analysis env)
- Test: `test_mecfs_bio/unit/build_system/wf/object_store/test_fake_object_store.py`

**Interfaces:**
- Produces: `ObjectHead(size_bytes:int, sha256:str|None)` frozen;
  `ObjectStore` ABC with `head(uri:str) -> ObjectHead | None` and `upload_from_url(source_url:str, uri:str) -> str` (streams URL→object with One Zone-IA + a stored SHA-256 checksum, returns the sha256);
  `S3ObjectStore(client=boto3.client("s3"))`; `FakeObjectStore(objects: dict[str, ObjectHead])` recording uploads in `.uploaded: list[tuple[str,str]]`.

- [ ] **Step 1: Write failing test for the fake** (`head` returns None for absent, the recorded `ObjectHead` for present; `upload_from_url` appends to `.uploaded` and makes the object subsequently `head`-able).
- [ ] **Step 2: Run → FAIL.**
- [ ] **Step 3: Implement** ABC + `FakeObjectStore` (in-memory), and `S3ObjectStore`: `head` via `head_object` (+ `get_object_attributes` with `ObjectAttributes=["Checksum"]` for the stored SHA-256, returning None if absent); `upload_from_url` streams via `upload_fileobj` with `ExtraArgs={"StorageClass":"ONEZONE_IA","ChecksumAlgorithm":"SHA256"}` and reads back the stored checksum. (S3ObjectStore has no unit test — it is exercised only through the fake and the manual first-staging run; document this.)
- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Commit.** `git commit -m "feat(objectstore): ObjectStore capability with S3 + fake"`

---

## Task 7: Wire capabilities into `WF`, `make_wf`, runner, config

**Files:**
- Modify: `mecfs_bio/build_system/wf/base_wf.py`
- Modify: `mecfs_bio/build_system/runner/simple_runner.py:115` (the `make_wf(...)` call)
- Modify: `mecfs_bio/analysis/runner/default_runner.py`
- Test: `test_mecfs_bio/unit/build_system/wf/test_make_wf_remote.py`

**Interfaces:**
- Consumes: `RemoteExecutor`, `SkyPilotRemoteExecutor`, `ObjectStore`, `S3ObjectStore`.
- Produces: `WF` gains `remote_executor: RemoteExecutor` and `object_store: ObjectStore`; `make_wf(..., remote_executor=None, object_store=None)` defaults them (SkyPilot + S3); `WF.remote_run(job, out_dir)` and `WF.object_store` accessor. New config keys read by `default_runner`: `remote_region`, `remote_s3_bucket`, `gctb_image`.

- [ ] **Step 1: Write failing test.**

```python
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.build_system.wf.remote_executor.fake_remote_executor import FakeRemoteExecutor

def test_make_wf_accepts_injected_remote_executor():
    fake = FakeRemoteExecutor()
    wf = make_wf(remote_executor=fake)
    assert wf.remote_executor is fake
```

- [ ] **Step 2: Run → FAIL.**
- [ ] **Step 3: Implement.** Add the two attributes to the frozen `WF` and params to `make_wf` (default `SkyPilotRemoteExecutor()` / `S3ObjectStore()`), add a `remote_run` convenience method. In `simple_runner.py` pass both through `make_wf(...)`. In `default_runner.py` add `_REMOTE_REGION_KEY`/`_REMOTE_BUCKET_KEY`/`_GCTB_IMAGE_KEY` readers (mirroring `_get_asset_root_path`) and construct the default `SkyPilotRemoteExecutor`/`S3ObjectStore` from them. Update the `default_runner.py` module docstring with the new optional config keys.
- [ ] **Step 4: Run → PASS**, and `pixi r invoke green` (import-linter: keep the new `wf` submodules within the existing layer).
- [ ] **Step 5: Commit.** `git commit -m "feat(wf): wire remote_executor + object_store capabilities"`

---

## Task 8: `StageGwfmReferenceTask`

**Files:**
- Create: `mecfs_bio/build_system/task/sbayesrc/stage_gwfm_reference_task.py`
- Create: `mecfs_bio/build_system/meta/gwfm_reference_marker_meta.py` (or reuse `SimpleFileMeta`; pick one and be consistent)
- Test: `test_mecfs_bio/unit/build_system/task/sbayesrc/test_stage_gwfm_reference_task.py`

**Interfaces:**
- Consumes: `GWFM_REFERENCE_BUNDLE`, `GWFM_REFERENCE_VERSION` (Task 1); `ObjectStore`, `ObjectHead` (Task 6); `FileAsset`, `Fetch`, `WF`.
- Produces: `StageGwfmReferenceTask(bucket: str)`; `execute` returns a `FileAsset` marker whose JSON is `{"version":..., "s3_prefix":"s3://<bucket>/sbayesrc/ld/<version>/", "files":[{"filename","size_bytes","sha256"}...]}` derived from the pinned bundle (deterministic across machines). Downstream reads `s3_prefix` + per-file remote layout.

- [ ] **Step 1: Write failing tests** (inject a `FakeObjectStore` via `wf`):
  - all bundle files already present+matching → `.uploaded == []`, marker written with the pinned manifest.
  - one file absent → that file's `upload_from_url` is called exactly once; others not re-uploaded.
  - marker JSON is byte-identical for two task instances with the same bucket (determinism).

```python
def test_skips_upload_when_bucket_populated(tmp_path, wf_with_fake_store_all_present):
    task = StageGwfmReferenceTask(bucket="mybucket")
    asset = task.execute(tmp_path, fetch=_noop_fetch, wf=wf_with_fake_store_all_present)
    assert wf_with_fake_store_all_present.object_store.uploaded == []
    assert '"version"' in Path(asset.path).read_text()
```

- [ ] **Step 2: Run → FAIL.**
- [ ] **Step 3: Implement.** For each `ReferenceBundleFile`: compute `uri = f"s3://{bucket}/sbayesrc/ld/{version}/{filename}"`; `head = wf.object_store.head(uri)`; if `head is None` or `head.size_bytes != f.size_bytes` or (`f.sha256` and `head.sha256 != f.sha256`) → `sha = wf.object_store.upload_from_url(f.source_url, uri)` (log the sha to record back into constants). Write the marker JSON (sorted keys) and return `FileAsset`. Assert the bundle is non-empty (fail fast).
- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Commit.** `git commit -m "feat(gwfm): StageGwfmReferenceTask with S3 dedup + deterministic marker"`

> After the FIRST real staging run, paste the recorded per-file SHA-256s into `GWFM_REFERENCE_BUNDLE` and commit, so future traces are checksum-verified.

---

## Task 9: `SumstatsToCojoMaTask`

**Files:**
- Create: `mecfs_bio/build_system/task/sbayesrc/sumstats_to_cojo_ma_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/sbayesrc/test_sumstats_to_cojo_ma_task.py`

**Interfaces:**
- Consumes: a sumstats parquet dep (read via `scan_dataframe_asset` per repo convention); `Fetch`, `WF`, `FileAsset`.
- Produces: `SumstatsToCojoMaTask(...)` writing a `.ma` with GCTB/COJO columns `SNP A1 A2 freq b se p N` (tab-separated, header row). Trait/project derived from the dep's meta (do NOT accept as args) via an isinstance assertion, per repo convention.

- [ ] **Step 1: Write failing test** on a tiny 3-row polars frame → assert the `.ma` has the 8 columns in order and correct values.
- [ ] **Step 2: Run → FAIL.**
- [ ] **Step 3: Implement** with polars: fetch dep, select/rename to the 8 COJO columns using column-name constants (add any missing to `constants/gwaslab_constants.py`), `write_csv(sep="\t")`. Return `FileAsset`.
- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Commit.** `git commit -m "feat(gwfm): SumstatsToCojoMaTask (parquet -> .ma)"`

---

## Task 10: `GctbFineMapTask`

**Files:**
- Create: `mecfs_bio/build_system/task/sbayesrc/gctb_fine_map_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/sbayesrc/test_gctb_fine_map_task.py`

**Interfaces:**
- Consumes: `SumstatsToCojoMaTask` (`.ma` FileAsset), `StageGwfmReferenceTask` (marker FileAsset); the CLI templates + resource constants (Task 1); `WF.remote_executor`; `DirectoryAsset`.
- Produces: `GctbFineMapTask(ma_task, reference_task, threads:int=DEFAULT_VCPUS)`; `deps` returns `[ma_task, reference_task]`; `execute` builds a `RemoteJob` (commands from the templates, `input_files={ma_path: PurePath("work/<trait>.ma")}`, `s3_inputs={marker.s3_prefix: PurePath("work/ref")}`, `output_files=[PurePath("work/out")]`, `resources=RemoteResources(DEFAULT_MEMORY_GB, threads, DEFAULT_DISK_GB, region)`), calls `wf.remote_executor.run(job, scratch_dir)`, and returns `DirectoryAsset(scratch_dir/"work/out")`.

- [ ] **Step 1: Write failing test** with `FakeRemoteExecutor` (stub the `out/` dir contents so `DirectoryAsset` post-init passes): assert the recorded `RemoteJob` has the 3 gctb commands in order, the `.ma` in `input_files`, the marker prefix in `s3_inputs`, and `resources.memory_gb == DEFAULT_MEMORY_GB`.

```python
def test_builds_wellformed_remote_job(tmp_path, ma_asset, marker_asset):
    fake = FakeRemoteExecutor(stub_outputs={PurePath("work/out/snpRes.txt"): "x"})
    wf = make_wf(remote_executor=fake)
    task = GctbFineMapTask(ma_task=_const(ma_asset), reference_task=_const(marker_asset))
    asset = task.execute(tmp_path, fetch=_fetch_map({...}), wf=wf)
    job = fake.last_job
    assert sum("gctb" in c for c in job.commands) == 3
    assert any(str(p).endswith(".ma") for p in job.input_files.values())
    assert job.resources.memory_gb == c.DEFAULT_MEMORY_GB
    assert Path(asset.path).is_dir()
```

- [ ] **Step 2: Run → FAIL.**
- [ ] **Step 3: Implement** (read the marker JSON via `fetch`, format the templates with concrete `work/`-relative paths, assemble the `RemoteJob`, call the executor, return the `DirectoryAsset`). If `USES_PRECOMPUTED_EIGEN`, drop the `--make-ldm-eigen` command.
- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Commit.** `git commit -m "feat(gwfm): GctbFineMapTask dispatches remote gctb job"`

---

## Task 11: GCTB Docker image + publish task

**Files:**
- Create: `docker/gctb/Dockerfile`
- Create: `docker/gctb/LICENSE-GCTB` (copy of GCTB's MIT license + copyright)
- Modify: `tasks.py` (add `build_push_gctb_image`)
- Test: `test_mecfs_bio/system/test_gctb_image.py`

**Interfaces:**
- Consumes: `GCTB_VERSION`, `GCTB_BINARY_URL`, `GCTB_BINARY_SHA256` (Task 1).
- Produces: a runnable public image tag `<registry>/gctb:<version>`; `invoke build-push-gctb-image` (maintainer-only).

- [ ] **Step 1: Write the Dockerfile.**

```dockerfile
FROM debian:stable-slim
RUN apt-get update && apt-get install -y --no-install-recommends libgomp1 ca-certificates curl unzip \
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

- [ ] **Step 2: Write the failing system test** (build the image locally, run `gctb` help in it, assert exit 0 and expected banner text present). Skip nothing — pixi guarantees Docker per repo convention.

```python
def test_gctb_image_runs(tmp_path):
    subprocess.run(["docker","build","--build-arg",f"GCTB_URL={c.GCTB_BINARY_URL}",
                    "--build-arg",f"GCTB_SHA256={c.GCTB_BINARY_SHA256}",
                    "-t","gctb:test","docker/gctb"], check=True)
    out = subprocess.run(["docker","run","--rm","gctb:test","gctb"], capture_output=True, text=True)
    assert "GCTB" in (out.stdout + out.stderr)
```

- [ ] **Step 3: Run → FAIL then implement** the `build_push_gctb_image` invoke task (builds with the two build-args, tags `<registry>/gctb:<GCTB_VERSION>`, pushes; registry from config). Document it as maintainer-only.
- [ ] **Step 4: Run → PASS.**
- [ ] **Step 5: Commit.** `git commit -m "feat(gwfm): minimal MIT-compliant GCTB image + publish task"`

---

## Task 12: End-to-end system test via `LocalDockerRemoteExecutor`

**Files:**
- Create: `test_mecfs_bio/system/test_gctb_gwfm.py`
- Create: `test_mecfs_bio/system/test_data/gwfm_toy/` (hello-world-scale LD + `.ma` + annotation + gene-map)

**Interfaces:**
- Consumes: `GctbFineMapTask`, `LocalDockerRemoteExecutor`, the real `gctb:test` image, toy data.

- [ ] **Step 1: Assemble toy inputs** — a tiny synthetic LD (a few hundred SNPs, one small block), matching `.ma`, annotation, gene-map, sized so `gctb --gwfm RC` finishes in seconds (mirror the existing "system test pattern for Docker-based GWAS tools"). Point `s3_inputs` at a local dir instead of S3 for this test (the LocalDockerRemoteExecutor treats an empty `s3_inputs` and mounts the toy reference via `input_files`).
- [ ] **Step 2: Write the test** — run `GctbFineMapTask.execute` with `make_wf(remote_executor=LocalDockerRemoteExecutor())`; assert the output dir contains the expected gctb result files and a non-empty PIP/credible-set file.
- [ ] **Step 3: Run → PASS** (`pixi r pytest test_mecfs_bio/system/test_gctb_gwfm.py -v`).
- [ ] **Step 4: Commit.** `git commit -m "test(gwfm): end-to-end GWFM via local docker executor on toy data"`

---

## Task 13: Collaborator docs + config template

**Files:**
- Modify: `default_runner_config.yaml` is gitignored — instead add a committed `default_runner_config.example.yaml` with the remote-exec keys commented.
- Create: a short setup section (where the repo keeps contributor docs — confirm location; the mkdocs `docs/Getting_Started/` tree is a candidate but is the published site, so a top-level `CONTRIBUTING`-style note or an `experiments/claude` note may be preferred).

- [ ] **Step 1** Write the example config (`remote_region`, `remote_s3_bucket`, `gctb_image`) with comments, plus the collaborator setup steps: `pixi install`, `aws configure`/SSO, `sky check`, fill config, run. Note `GWFM_ASSUME_YES=1` for non-interactive runs.
- [ ] **Step 2** `pixi r invoke green`, then commit. `git commit -m "docs(gwfm): collaborator setup + example remote config"`

---

## Self-Review

**Spec coverage:** Approach-1 WF seam → Tasks 2–7,10; SkyPilot/AWS → Task 5; public MIT image → Task 11; reference bundle (LD+eigen+annot+gene-map) staging with deterministic marker + S3 checksum dedup → Task 8; `.ma` conversion → Task 9; on-demand-only / no `use_spot` → Task 2 note + Global Constraints; cost/failure safety (teardown finally, idle-autostop, pre-launch confirm, streamed logs) → Task 5; One Zone-IA same-region → Task 6/Constraints; testing (fake unit + local-docker system) → Tasks 3,4,10,12; config/credentials (AWS chain, gitignored config) → Tasks 7,13; future local-annotation asset → out of scope (recorded in spec). Recon-dependent externals → Task 1.

**Placeholders:** the only deferred values (bundle SHA-256s, exact gctb flags, `USES_PRECOMPUTED_EIGEN`, disk size, SkyPilot output-retrieval call) are explicitly produced by Task 1 / recorded after first staging, not left vague in consumer tasks.

**Type consistency:** `RemoteJob`/`RemoteResources` fields, `RemoteExecutor.run(job, local_output_dir)`, `ObjectStore.head/upload_from_url`, `ObjectHead(size_bytes, sha256)`, and the marker JSON schema (`version`/`s3_prefix`/`files`) are used identically across Tasks 2–12.

## Open verification items carried from the spec (resolved in Task 1)
1. GCTB license → **MIT, redistribution OK** (already resolved; image includes LICENSE-GCTB).
2. GCTB resume support → set `GCTB_SUPPORTS_RESUME`; gates any future spot work.
3. Precomputed `eigen/` usable directly → `USES_PRECOMPUTED_EIGEN`; VM disk sizing.
4. Exact GWFM annotation + gene-map files/sizes on gctbhub.
5. Current S3 One Zone-IA pricing for the chosen region.
6. Exact SkyPilot output-retrieval API for the installed version (Task 5 Step 3).
