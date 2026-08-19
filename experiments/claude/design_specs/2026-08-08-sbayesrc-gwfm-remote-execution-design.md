# SBayesRC Genome-Wide Fine-Mapping via Remote (Cloud) Execution — Design

Date: 2026-08-08
Status: Approved design, pending implementation plan

## Goal

Add a build-system task that runs SBayesRC **genome-wide fine-mapping (GWFM)** on a
trait's GWAS summary statistics. A single global run over ~13 million SNPs needs
**~150 GB RAM, ~24 CPU cores, ~13 h** (per Wu et al., *Nat Genet* 2026), far beyond the
16 GB development laptop. The task must therefore run the heavy step on a transient
cloud instance while the rest of the build graph runs locally, preserving the build
system's caching and lineage guarantees.

Usage profile: a **handful of traits**, one global run each (GWFM is not per-locus),
run infrequently.

## Key facts established during brainstorming

- **GWFM is pure `gctb` binary** — no R, no rpy2, no SBayesRC R package. Three CLI calls:
  1. `gctb --make-ldm-eigen` — match the LD reference to the trait's sumstats.
  2. `gctb --gwfm RC ...` — the heavy ~13 h / ~150 GB MCMC.
  3. `gctb --cs ...` — credible sets (cheap).
- **No suitable off-the-shelf GCTB image exists.** `jianzeng/GCTB` has no Dockerfile/image
  and no in-repo binary (compile-from-source: C++11, OpenMP, Eigen3, Boost). The gctbhub
  software page ships a precompiled `gctb_2.5.5_Linux.zip` (GWFM needs `gctb >= 2.5.1`).
  The existing `zhiliz/sbayesrc` image tops out at `0.2.6` (the SBayesRC R package for
  polygenic prediction) and predates GWFM. Third-party GCTB images are unofficial/unpinned.
- **LD reference (UKB EUR Imputed 13M)** from gctbhub is **~200 GB**, immutable, and identical
  for every run. Accessed ~1–2×/month.

## Decisions

1. **Architecture — remote *command* executor as a `WF` capability** (not full
   build-system-managed remote execution, not an out-of-band script). The task's Python
   `execute` runs locally for light orchestration; only the `gctb` commands run in the
   cloud container. Rationale: the heavy unit of work is already a containerized command,
   so shipping just that command (rather than reproducing the whole pixi/R/Docker
   environment in the cloud) is far simpler and keeps caching/lineage intact. The task is
   agnostic to the provisioning *mechanism* (it only builds a `RemoteJob`), but explicitly
   opts its heavy step into remote execution.
2. **Provider + provisioning:** AWS via **SkyPilot**. Standard AWS credential chain (no
   bespoke credentials file).
3. **Container image:** a minimal **self-built public image** (`debian-slim` + `libgomp1` +
   the pinned, checksummed `gctb` binary). Committed `Dockerfile` + a maintainer-only
   `invoke build-push-gctb-image`. Public so collaborators pull with zero auth. Pending a
   check of GCTB's redistribution terms; if unclear, the Dockerfile downloads the
   checksummed binary from gctbhub at build time and each user builds/pushes to their own
   registry.
4. **Instance type:** on-demand only. Spot is deferred and will be considered **only** if
   GCTB's docs/code show explicit MCMC resumption support. (`use_spot` is therefore *not*
   part of the interface yet.)
5. **Reference data storage:** shared **S3 One Zone-IA** bucket in the **same region as
   compute** (free S3→EC2 transfer; single-AZ durability acceptable because the data is
   re-derivable from gctbhub). ~$4–6/mo.
6. **Cost/failure safety:** guaranteed teardown, idle-autostop backstop, pre-launch
   confirmation with cost estimate, streamed remote logs (details below).
7. **Testing:** fake executor for unit tests; a `LocalDockerRemoteExecutor` running the
   identical container locally on toy data for a single system test. No 13 h / ~$50 cloud
   run ever hits CI.

## Data flow

```
LOCAL (build graph)                                  CLOUD (ephemeral, per-run)
─────────────────────────────                        ──────────────────────────
trait sumstats (parquet)
   │  SumstatsToCojoMaTask  (light, local)
   ▼
 trait.ma  ───────────────────────ship (small)───►  /work/trait.ma
                                                          │  VM setup: pull LD ref
S3 bundle marker ◄─ StageGwfmReferenceTask                │  from S3, docker-pull image
   (LD + eigen + annot + gene-map,                         │
    points at s3://…/Imputed13M/v1/, checksummed)          │
                                                          ▼
GctbFineMapTask.execute():                          docker run gctb:
  build 3 gctb commands ─► wf.remote_executor.run() ──►  1) --make-ldm-eigen
                                                         2) --gwfm RC  (~13 h)
  local scratch_dir ◄──── retrieve outputs ◄──────────   3) --cs
       │                   (finally: sky down)
       ▼
  FineMapResult asset ─► rebuilder hashes/caches normally
```

Because `GctbFineMapTask` returns a normal **local** asset, the verifying-trace rebuilder
hashes and caches it as usual: an unchanged trait launches no instance.

## Component 1 — the remote-executor seam

New `WF` capability, alongside `downloader` / `synapse_downloader`. Typed `@frozen` job
description (per repo convention of named-attribute objects over bare tuples):

```python
@frozen
class RemoteResources:
    memory_gb: int
    vcpus: int
    disk_gb: int
    region: str | None = None          # from gitignored config

@frozen
class RemoteJob:
    image: str                              # digest-pinned public gctb image
    commands: Sequence[str]                 # gctb calls, run in-container on the VM
    input_files: Mapping[Path, PurePath]    # local -> remote (small per-run inputs)
    s3_inputs: Mapping[str, PurePath]       # s3 uri -> remote (large shared refs, pulled on VM)
    output_files: Sequence[PurePath]        # remote paths to retrieve
    resources: RemoteResources

# Explicit named-instance override (e.g. r6i.8xlarge) is deliberately NOT supported
# yet. If added later, refactor RemoteJob.resources to `RemoteResources | ExplicitInstance`
# so that forcing a named type does not redundantly require vcpus/memory_gb — keeping
# invalid states unrepresentable rather than adding an optional `instance_type` field.

class RemoteExecutor(Protocol):
    def run(self, job: RemoteJob, local_output_dir: Path) -> None: ...
```

- **Production:** `SkyPilotRemoteExecutor` — builds a SkyPilot task (`file_mounts` for small
  inputs; `setup:` pulls the S3 reference + `docker pull`s the image and verifies the
  reference exists, failing fast if missing; `run:` executes `commands` in the container),
  launches on-demand, syncs `output_files` back into `local_output_dir`, then tears down.
- **System-test:** `LocalDockerRemoteExecutor` — runs the identical container locally,
  ignoring `resources`, so the real `gctb` is exercised end-to-end on toy data.
- **Unit-test:** a fake executor that records the `RemoteJob` for assertions.

`make_wf` gains a `remote_executor` parameter defaulting to `SkyPilotRemoteExecutor`.

## Component 2 — `StageGwfmReferenceTask`

Stages the full immutable **GWFM reference bundle** to the shared S3 bucket; the **local
asset is a small deterministic marker**, so we never pull the bundle locally.

The bundle is everything GWFM needs from gctbhub that is shared and immutable across all
traits:
- the ~200 GB UKB EUR Imputed 13M LD matrix (+ `eigen/`),
- the **SNP annotation** matched to the 13M reference (`--annot`), and
- the **gene map** (`--gene-map`).

These are staged and pulled together because they are all shared/immutable gctbhub
resources; folding the annotation (potentially multi-GB) into the bundle avoids shipping it
from the laptop on every run. Note: the repo's existing `baseline_2_2_annotation` asset is
HapMap3-scale for the polygenic-prediction workflow and is **not** the 13M GWFM annotation —
the correct GWFM annotation file must be identified on gctbhub (see open items).

**Future, separate asset (not part of this bundle):** a small task that fetches a *local*
copy of the annotation for building plots/tables that illustrate how annotations shaped the
fine-mapping results. This is deliberately decoupled from `StageGwfmReferenceTask` so
plotting never forces the 200 GB bundle onto the laptop; it can share the same pinned
manifest/version to stay consistent.

- **Pinned in code:** a manifest of `(filename, size, sha256)` for every file in the bundle
  (LD, eigen, annotation, gene map), plus a `version` string (e.g. `Imputed13M/v1`). S3
  target key includes the version and is never overwritten.
- **`execute()`:**
  1. For each manifest file, check S3 **without downloading**: `HEAD` (size) +
     `GetObjectAttributes` (stored checksum).
  2. All present + matching → no S3 work.
  3. Missing/mismatched → stream from gctbhub → upload to S3 with the checksum recorded.
     (Only the first stager pays this.)
  4. Write the local marker JSON (S3 prefix + `version` + manifest); return it as the
     `FileAsset`.
- **Trace stability:** marker content derives from the *pinned manifest + configured S3
  location*, not live S3 state, so it hashes identically across machines → stable downstream
  dependency trace → good caching. Bumping `version` on an upstream republish forces a
  correct rebuild.
- **Collaborator semantics:** a fresh clone has no local marker → `execute()` runs → a few
  `HEAD`/`GetObjectAttributes` calls confirm the shared bucket is already populated → **no
  re-upload**, just writes the local marker (seconds).

## Component 3 — supporting tasks

| Task | Runs | Notes |
|---|---|---|
| `SumstatsToCojoMaTask` | local | trait parquet → gctb `.ma`; near-existing (`build_protein_cojo_ma.py` on the `polypwas_s_bayes_tasks` branch) |
| `StageGwfmReferenceTask` | local, idempotent | Component 2 — LD + eigen + annotation + gene map |
| `GctbFineMapTask` | local orchestration → remote command | heavy step; deps = `.ma`, GWFM reference-bundle marker |
| `build-push-gctb-image` | `invoke` task (not a build Task) | maintainer-only infra |

## Cost / failure safety

- **Guaranteed teardown:** the executor wraps launch in `try/finally` so `sky down` runs even
  on exception / Ctrl-C.
- **Idle-autostop backstop:** SkyPilot idle-autostop reaps the VM if the local process dies.
- **Pre-launch confirmation:** print instance type + cost estimate and require confirmation
  before provisioning; overridable by a config flag / env var for non-interactive runs.
- **Streamed logs:** the remote job log streams locally so progress is visible across ~13 h.

## Collaborator setup (target: minimal)

`pixi install` (gets SkyPilot + Docker client) → `aws configure` (or SSO) once → `sky check`
→ put region + S3 bucket into the gitignored `default_runner_config.yaml` → run the task.
No secrets in the repo; no image to build (they consume the public one).

## Configuration

Non-secret, user-specific settings live in the existing gitignored `default_runner_config.yaml`
(where `path_remap` already lives): AWS region, S3 reference-bucket name, default instance
type, and the digest-pinned image reference. Sensible defaults committed in code; the local
file overrides only what's user-specific. AWS credentials come from the standard AWS chain,
never from the repo.

## Testing strategy

- **Unit:** inject the fake executor; assert the `RemoteJob` is well-formed (commands, mounts,
  declared outputs). No cloud, no Docker.
- **System:** `LocalDockerRemoteExecutor` runs the real `gctb` on hello-world-scale toy LD /
  sumstats, mirroring the existing system-test pattern for Docker-based GWAS tools.

## Open items to resolve during implementation

1. Confirm GCTB's redistribution terms → decide public prebuilt image vs. build-time download.
2. Check GCTB docs/code for MCMC resumption → gate any future spot support.
3. Confirm whether GWFM consumes the precomputed `eigen/` directory directly (skipping
   `--make-ldm-eigen`) and finalize the exact reference file manifest + VM disk sizing
   (~500 GB `gp3` planned).
4. Identify the exact GWFM SNP-annotation file and gene-map on gctbhub (matched to the 13M
   reference; distinct from the repo's HapMap3 `baseline_2_2_annotation`), and their sizes.
5. Confirm current S3 One Zone-IA pricing for the chosen region.
```
