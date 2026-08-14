# Running remote GWFM (SBayesRC / GCTB genome-wide fine-mapping)

These notes let a collaborator run the GCTB genome-wide fine-mapping (GWFM) tasks,
which dispatch a container to a transient on-demand AWS instance via SkyPilot and
retrieve the result. Everything runs from the local build system; the AWS instance
exists only for the duration of a job and is torn down afterwards.

## One-time setup

1. **Install the environment.** From the repo root:

   ```
   pixi install
   ```

   This provides pixi-managed tools including the SkyPilot and AWS CLIs; run
   everything through `pixi r` so it resolves to the pinned versions.

2. **Configure AWS credentials.** Credentials are read only from the standard AWS
   credential chain (environment variables, `~/.aws/`, SSO, or an instance role) —
   never from anything checked into the repo. Set up whichever your organization
   uses, for example:

   ```
   pixi r aws configure          # long-lived access key + region
   # or, for SSO:
   pixi r aws configure sso
   pixi r aws sso login
   ```

   Confirm they work: `pixi r aws sts get-caller-identity`.

3. **Verify SkyPilot can see AWS.**

   ```
   pixi r sky check aws
   ```

   This must report AWS as enabled before any remote job will launch.

The AWS region, the shared S3 reference bucket, and the GCTB container image are
properties of the tasks themselves, not per-user settings: the reference bucket and
image are shared across all collaborators, and compute runs in the region of the
reference data (where the bulk of the data transfer is), so an individual user does
not choose them. The one remote-exec setting a collaborator must supply is the
scratch prefix, configured in the machine-local runner config (below).

## Runner config for a run

Remote-exec settings live in the gitignored `default_runner_config.yaml` (copy it
from `default_runner_config.example.yaml`), alongside the general asset-store keys
(`asset_root`, `info_store`, `path_remap`). The default runner reads them and builds
the SkyPilot executor accordingly:

- `remote_scratch_s3` — **required for remote runs.** An `s3://` prefix the executor
  uses for scratch (job outputs staged through S3), e.g.
  `remote_scratch_s3: s3://my-gwfm-bucket/remote-exec-scratch`. There is no default;
  a remote task started without it fails fast at launch with a message naming this
  key.
- `remote_non_interactive` — **optional, defaults to false.** When true, the executor
  auto-approves each paid launch instead of prompting. Set it only in unattended
  contexts (CI, batch runs); leave it false so an interactive run always reviews the
  instance type and estimated cost before spending money.

## Running

With AWS reachable and `remote_scratch_s3` set, run the GWFM task the same way as any
other build-system task (through `pixi r`). Before launch the executor prints the
chosen instance type and an estimated cost; unless `remote_non_interactive: true` is
set it waits for confirmation. The instance is torn down when the job finishes,
including on failure.

## Cost & safety notes

- Instances are on-demand only (no spot), sized by the task's `RemoteResources`.
- The executor sets an idle-autostop so a stuck or abandoned instance does not run
  indefinitely, and tears the instance down in a `finally` block.
- Reference staging uses S3 with checksum-based dedup, so re-running a job does not
  re-upload the (large) reference bundle.
