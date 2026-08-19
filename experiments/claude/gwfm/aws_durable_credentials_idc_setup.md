# Durable local AWS credentials via IAM Identity Center

## Why this exists

Staging the ~192 GB GWFM LD reference to S3 is a ~10 h upload, and the later SBayesRC
genome-wide fine-mapping run is supervised locally for ~13 h. Both need local AWS
credentials that stay valid for the whole time.

The `aws login` flow used previously is **hard-capped at 12 hours per login** (AWS
refreshes the cached credentials every 15 min but the overall session cannot exceed the
IAM principal's max session duration, 12 h), after which it fails with
`LoginRefreshRequired` / `ExpiredToken` and needs a manual re-login. That killed the
upload mid-transfer and cannot cover the 13 h fine-mapping run at all.

IAM Identity Center (IdC) fixes this at the root: the CLI obtains temporary credentials
that **auto-refresh silently for up to the configured access-portal session duration**
(15 min – 90 days), with no browser re-auth in between. Set that duration above the
longest run and the credential-expiry problem disappears for both the upload and the
SBayesRC supervision — with no build-system code changes. It also moves day-to-day work
off the account root user, which is the recommended posture.

### How the refresh works (two separate timers)

- **Permission-set session** (1–12 h): the actual STS role credentials. The SDK/CLI
  refreshes these transparently whenever they near expiry, so 12 h here is not a wall —
  it is just the interval between silent refreshes. Set it to the 12 h maximum.
- **Access-portal (IdC) session** (15 min – 90 days, default 8 h): how long before the
  browser sign-in must be repeated. With the refresh-token OIDC client that
  `aws configure sso` sets up (`sso_registration_scopes = sso:account:access`), the CLI
  silently mints new access tokens across this whole window. **This is the timer that
  must exceed your longest unattended run** — the default 8 h does NOT cover the 13 h
  SBayesRC job, so raise it.

## Prerequisites

- AWS CLI v2 ≥ 2.32.0 (`pixi r aws --version`; PKCE SSO needs ≥ 2.22.0). Confirmed
  present at 2.36.x.
- Admin (or root) access to the AWS account to enable IdC and create the user /
  permission set. This is a one-time console step.

## Part 1 — Console, one-time (as account admin/root)

Enabling IdC on a standalone account automatically creates an AWS Organization
containing just this account (free), and pins an IdC "home Region" — choose the Region
you already use (`us-east-1`); it is effectively permanent.

1. **Enable IAM Identity Center.** Console → *IAM Identity Center* → *Enable*. Pick the
   home Region (`us-east-1`).

2. **Create a user.** IdC → *Users* → *Add user* (e.g. username `gwfm-ops`, your email).
   Accept the email invitation and set a password (this is your new non-root identity).

3. **Create a permission set.** IdC → *Permission sets* → *Create*.
   - Predefined `AdministratorAccess` is the simplest choice and is strictly better than
     the root usage it replaces. (Least-privilege alternative below.)
   - Set the **permission-set session duration to 12 hours** (the maximum).

4. **Assign the user to the account.** IdC → *AWS accounts* → select this account →
   *Assign users or groups* → pick `gwfm-ops` → attach the permission set from step 3.

5. **Raise the access-portal session duration.** IdC → *Settings* → *Authentication* /
   *Session settings* → set the session duration **longer than your longest run**.
   - Minimum useful value here is > 13 h (to cover the SBayesRC run). A value like
     **24 h**, or up to **7 days** for convenience across a multi-day
     staging-then-run workflow, is reasonable. Longer = fewer re-logins but a longer
     window before forced re-auth; lower it later if you prefer. Max is 90 days.

### Least-privilege alternative to AdministratorAccess

SkyPilot's *first ever* launch creates the shared `skypilot-v1` IAM role + instance
profile, which needs IAM write actions (`iam:CreateRole`, `CreateInstanceProfile`,
`AttachRolePolicy`, `AddRoleToInstanceProfile`, `PassRole`). `PowerUserAccess` excludes
IAM, so with it you must either (a) pre-create the `skypilot-v1` instance profile once as
admin, after which PowerUserAccess suffices for all later launches, or (b) attach a small
custom policy granting exactly those IAM actions. For a single-maintainer research
account, `AdministratorAccess` avoids this friction.

## Part 2 — Laptop, one-time (replace the old login profile)

Mixed credential types in one profile cause `ExpiredToken` errors, so clear the old
`aws login` session first.

1. **Sign out the old session:**

   ```
   pixi r aws logout
   ```

   Also remove the old `[default]` block containing `login_session = ...` from
   `~/.aws/config` (the wizard below will rewrite `[default]`).

2. **Configure the SSO profile as the default:**

   ```
   pixi r aws configure sso
   ```

   Answer the prompts:
   - `SSO session name`: e.g. `gwfm`
   - `SSO start URL`: from IdC → *Settings* → *Summary* (or *Dashboard*), the AWS access
     portal URL (looks like `https://d-xxxxxxxxxx.awsapps.com/start`)
   - `SSO region`: `us-east-1` (the IdC home Region)
   - `SSO registration scopes`: `sso:account:access`
   - Select the account and the permission set (role) when listed.
   - `Default client Region`: `us-east-1`
   - `Profile name`: **`default`** (so SkyPilot, boto3, and the staging script all pick
     it up with no `AWS_PROFILE` needed)

   The resulting `~/.aws/config` should look like:

   ```
   [profile default]
   sso_session = gwfm
   sso_account_id = <ACCOUNT_ID>
   sso_role_name = AdministratorAccess
   region = us-east-1
   output = json

   [sso-session gwfm]
   sso_region = us-east-1
   sso_start_url = https://d-xxxxxxxxxx.awsapps.com/start
   sso_registration_scopes = sso:account:access
   ```

3. **Log in and verify durability:**

   ```
   pixi r aws sso login
   pixi r aws configure list          # TYPE column must read `sso`, not `login`
   pixi r aws sts get-caller-identity  # shows the assumed-role ARN, not :root
   ```

## Part 3 — SkyPilot + repo

- Confirm SkyPilot sees the credentials:

  ```
  pixi r sky check aws
  ```

- SkyPilot does **not** upload SSO credentials to instances; it attaches the
  `skypilot-v1` IAM instance profile, whose role credentials come from the instance
  metadata service and auto-refresh indefinitely (they never expire). So on-instance
  S3 work during a job is unaffected by any local session limit. With
  `AdministratorAccess`, the first-launch instance-profile creation just works.
- The repo's `default_runner_config.yaml` region is already `us-east-1`; no change
  needed. Credentials are read only from the standard AWS chain, never from the repo.

## Part 4 — Re-run the staging upload

1. Abort any leftover incomplete multipart upload first (see
   `stage_gwfm_reference_s3_setup.md`); the 7-day `AbortIncompleteMultipartUpload`
   lifecycle rule is a backstop but clearing it now stops the storage charge.
2. With the access-portal session duration set above ~11 h, re-run staging; the
   committed multipart-chunk fix keeps the 192 GB object within the 10,000-part limit,
   and durable credentials let the ~10 h transfer complete without a mid-run expiry:

   ```
   pixi r python experiments/tralfamadorian97/gwfm/stage_gwfm_ref.py
   ```

   Staging restarts the big file from zero (no resume); already-present bundle files are
   skipped by the size/checksum dedup check.

## Re-authenticating later

When the access-portal session finally lapses, a single `pixi r aws sso login` restores
durable credentials; there is no per-command or 12 h interruption in between. Use
`pixi r aws logout` to end the session deliberately.
