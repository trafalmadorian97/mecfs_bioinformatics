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

2. **Configure durable AWS credentials.** Credentials are read only from the standard
   AWS credential chain (environment variables, `~/.aws/`, SSO, or an instance role) —
   never from anything checked into the repo.

   **Use a credential that does not expire mid-run.** The GWFM fine-mapping job is
   supervised locally for ~13 h, so the browser `aws login` flow will **not** work — it
   is hard-capped at 12 h per login and dies mid-run with `ExpiredToken`. (Staging the
   reference bundle is a one-time owner action; as a collaborator you only ever do quick
   `HEAD` checks against the shared bucket, so the credential-lifetime concern is really
   about the long fine-mapping run.) Pick whichever fits your account:

   | Your situation | What to do | Notes |
   |---|---|---|
   | Your org already uses AWS SSO / IAM Identity Center | `pixi r aws configure sso`, then `pixi r aws sso login` | Simplest — your admin already set up IdC; just use your org's start URL. |
   | A personal / standalone AWS account | Static IAM user access keys (detailed below) | Never expire; no IdC to enable. |
   | You were given an IdC user in the data owner's account | `pixi r aws configure sso` with the owner's start URL, then `pixi r aws sso login` | The owner does the one-time user + permission-set assignment. |

   Whichever you choose, the identity needs, in the account you launch into: **EC2 +
   IAM** (so SkyPilot can create/use the shared `skypilot-v1` instance profile on the
   first launch) and **S3 read**. `AdministratorAccess` covers this; `PowerUserAccess`
   works once the `skypilot-v1` instance profile already exists.

   Confirm whatever you set up works:

   ```
   pixi r aws sts get-caller-identity
   ```

   ### Static IAM user access keys (simplest for a personal account)

   Long-lived access keys never expire, so they sidestep the session-length problem
   entirely. The trade-off is a long-lived secret stored on disk — scope it, and delete
   it when you are done. In the AWS account you will launch into:

   1. **Create an IAM user.** Console → *IAM* → *Users* → *Create user*, e.g.
      `gwfm-runner`. Leave *console access* off (programmatic use only).
   2. **Attach permissions.** On *Set permissions* choose *Attach policies directly* and
      attach `AdministratorAccess` (simplest), or `PowerUserAccess` plus the IAM actions
      SkyPilot needs on its first launch: `iam:CreateRole`, `iam:CreateInstanceProfile`,
      `iam:AttachRolePolicy`, `iam:AddRoleToInstanceProfile`, `iam:PassRole`,
      `iam:GetRole`, `iam:GetInstanceProfile`.
   3. **Create an access key.** Open the user → *Security credentials* → *Create access
      key* → *Command Line Interface (CLI)* → create. Copy the **Access key ID** and
      **Secret access key** now — the secret is shown only once.
   4. **Configure the CLI locally:**

      ```
      pixi r aws configure
      ```

      Paste the Access key ID and Secret access key, set the default region to the
      reference-data region (`us-east-1`), and output `json`. This writes the keys to
      `~/.aws/credentials` and the region to `~/.aws/config`.
   5. **Verify:**

      ```
      pixi r aws sts get-caller-identity   # shows the gwfm-runner user ARN
      ```

   To rotate or revoke later: IAM → the user → *Security credentials* → deactivate or
   delete the access key (and delete the user once the project is done).

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
