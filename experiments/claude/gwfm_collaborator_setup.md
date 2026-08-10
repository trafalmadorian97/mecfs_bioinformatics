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

4. **Fill in the machine-local config.** Copy the committed example to the real
   (gitignored) config and edit it:

   ```
   cp default_runner_config.example.yaml default_runner_config.yaml
   ```

   Set the three remote-exec keys:

   - `remote_region` — the AWS region for compute and the S3 reference bucket. Keep
     the bucket and compute in the same region so reference transfer stays
     intra-region.
   - `remote_s3_bucket` — the bucket holding the staged GWFM reference bundle and
     scratch data.
   - `gctb_image` — the public GCTB container image reference. Build and push it
     once with `pixi r invoke build-push-gctb-image --registry <your-registry>`
     (see `docker/gctb/`), then point this key at the resulting
     `<your-registry>/gctb:<version>` tag.

   `default_runner_config.yaml` stays gitignored — do not commit it or put
   credentials in it.

## Environment variables for a run

The SkyPilot executor reads two environment variables at run time:

- `REMOTE_EXEC_SCRATCH_S3` — **required.** An `s3://` prefix the executor uses for
  scratch (job inputs/outputs staged through S3), e.g.
  `REMOTE_EXEC_SCRATCH_S3=s3://my-gwfm-bucket/remote-exec-scratch`.
- `REMOTE_EXEC_ASSUME_YES=1` — **optional.** Skips the interactive pre-launch cost
  confirmation prompt. Set it in non-interactive contexts (CI, batch runs); leave
  it unset to review the instance type and estimated cost before each launch.

## Running

With the config filled in and AWS reachable, run the GWFM task the same way as any
other build-system task (through `pixi r`). Before launch the executor prints the
chosen instance type and an estimated cost; unless `REMOTE_EXEC_ASSUME_YES=1` is
set it waits for confirmation. The instance is torn down when the job finishes,
including on failure.

## Cost & safety notes

- Instances are on-demand only (no spot), sized by the task's `RemoteResources`.
- The executor sets an idle-autostop so a stuck or abandoned instance does not run
  indefinitely, and tears the instance down in a `finally` block.
- Reference staging uses S3 with checksum-based dedup, so re-running a job does not
  re-upload the (large) reference bundle.
