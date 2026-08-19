# Creating the S3 scratch bucket for remote GWFM runs

Step-by-step setup for the S3 **scratch** bucket the SkyPilot remote executor uses to
shuttle a remote job's outputs back to your machine. This is a companion to
`stage_gwfm_reference_s3_setup.md`, but the two buckets are deliberately different
animals — do **not** copy the reference bucket's access model here (see "Why not reuse
the reference bucket" below).

## What the scratch bucket is for

`SkyPilotRemoteExecutor` (in
`mecfs_bio/build_system/wf/remote_executor/skypilot_remote_executor.py`) has no clean
SDK API for recursively downloading a finished job's outputs, so it round-trips them
through S3:

1. `run()` derives a run-scoped prefix `output_s3_prefix = {scratch_s3}/{cluster}`,
   where `scratch_s3` is the executor's configured scratch prefix and `cluster` is a
   fresh per-run id.
2. On the remote instance, the job's `run` phase uploads each output there with
   `aws s3 cp --recursive <output> <output_s3_prefix>/...` (see `build_sky_task`).
3. Back on your machine, `_retrieve_outputs` mirrors that prefix down with
   `aws s3 cp --recursive <output_s3_prefix>/... <local_dir>`.

`scratch_s3` comes from the `remote_scratch_s3` key in the machine-local
`default_runner_config.yaml`. Without it, `run()` fails fast at launch with
"scratch_s3 is not set".

## Why not reuse the reference bucket (or its setup)

The reference bucket (`stage_gwfm_reference_s3_setup.md`) is **public-read + Requester
Pays** so any collaborator can pull the 192 GiB LD reference and pay their own egress.
The scratch bucket is the opposite: it holds *your own* transient run outputs, read and
written only by you and your own compute.

The one setting that would actively break a run: **Requester Pays**. Only the
reference *read* passes `--request-payer requester` (the reference bucket needs it).
The scratch **write** (on the instance) and the scratch **read** (`_retrieve_outputs`
on your laptop) do **not** pass that flag, so a Requester-Pays scratch bucket would
reject both with 403. Keep the scratch bucket a plain private bucket with Payer =
BucketOwner.

| Aspect | Reference bucket | Scratch bucket |
|---|---|---|
| Region | compute region (`us-east-1`) | compute region (`us-east-1`) |
| Requester Pays | yes | **no** (would 403 the run) |
| Public read policy (`Principal:"*"`) | yes | **no** — fully private |
| Block Public Access | relaxed for the `"*"` policy | all four flags on |
| Lifecycle | abort-incomplete-MPU | abort-incomplete-MPU **+ expire old scratch** |
| Extra IAM policy | least-priv staging policy (optional) | none needed |

## Step 1: Region

Put the scratch bucket in the **same region as your GWFM compute** (`us-east-1`), the
same region as the reference bucket. The scratch round-trip is instance -> S3 -> laptop;
co-locating bucket and compute keeps the instance -> S3 leg intra-region. (Fine-mapping
outputs are small, so this matters far less than for the reference, but there is no
reason to split regions.)

## Step 2: Create the bucket (private)

Pick a globally-unique, DNS-compliant name. This repo uses
`mecfs-bio-remote-exec-scratch`.

```
pixi r aws s3api create-bucket --bucket mecfs-bio-remote-exec-scratch --region us-east-1
```

That is all the access configuration needed. Leave **Block Public Access fully on**
(the account/bucket default) and do **not** enable Requester Pays or attach any bucket
policy — the bucket should stay private.

## Step 3: Lifecycle (auto-clean transient scratch)

Scratch accumulates: one prefix per run, plus any parts orphaned by a hard kill. Two
rules keep it from costing anything over time — abort incomplete multipart uploads
after 7 days, and expire finished scratch objects after 30 days. Save as
`scratch-lifecycle.json`:

```json
{
  "Rules": [
    {
      "ID": "AbortIncompleteMPU",
      "Status": "Enabled",
      "Filter": {"Prefix": ""},
      "AbortIncompleteMultipartUpload": {"DaysAfterInitiation": 7}
    },
    {
      "ID": "ExpireOldScratch",
      "Status": "Enabled",
      "Filter": {"Prefix": ""},
      "Expiration": {"Days": 30}
    }
  ]
}
```

```
pixi r aws s3api put-bucket-lifecycle-configuration \
  --bucket mecfs-bio-remote-exec-scratch \
  --lifecycle-configuration file://scratch-lifecycle.json
```

If your runs keep outputs you want to hold longer than 30 days, raise or drop the
`ExpireOldScratch` rule — but remember these are meant to be transient copies; the
authoritative outputs land in your local asset store after `_retrieve_outputs`.

## Step 4: Permissions (nothing to attach)

No extra IAM policy is required, as long as the bucket, the compute, and you are all in
the same AWS account:

- The remote instance reads/writes scratch under the SkyPilot-managed `skypilot-v1`
  instance role, which carries `AmazonS3FullAccess`.
- Your local machine reads scratch under whatever identity the standard AWS credential
  chain resolves (here the `mecfs-bio-static-runner` IAM user, which has broad S3
  access).

Both principals live in account `920441304540`, so both already have full access to a
bucket in that account. (A collaborator running from a locked-down identity would
instead need `s3:GetObject`/`s3:PutObject`/`s3:ListBucket` on their own scratch bucket;
the owner does not.)

## Step 5: Wire it into the runner config

Add the prefix to your gitignored `default_runner_config.yaml` (copy the commented
example from `default_runner_config.example.yaml`):

```yaml
remote_scratch_s3: s3://mecfs-bio-remote-exec-scratch/remote-exec-scratch
```

A trailing sub-prefix (`/remote-exec-scratch`) is optional but tidy — every run then
nests under it as `.../remote-exec-scratch/remote-exec-<id>/`.

Leave `remote_non_interactive` unset (defaults to false) so an interactive run still
prints the instance type and estimated cost and waits for your confirmation before
spending money.

## Step 6: Verify

```
# Lifecycle rules present and Enabled
pixi r aws s3api get-bucket-lifecycle-configuration --bucket mecfs-bio-remote-exec-scratch
# Block Public Access: all four flags true
pixi r aws s3api get-public-access-block --bucket mecfs-bio-remote-exec-scratch
# Payer must be BucketOwner (NOT Requester)
pixi r aws s3api get-bucket-request-payment --bucket mecfs-bio-remote-exec-scratch
```

## Cost notes

- **Storage:** scratch holds small fine-mapping outputs briefly, and the 30-day
  expiration rule reaps them, so ongoing storage is negligible.
- **Transfer:** instance -> S3 (upload) is free; S3 -> laptop (download) is normal S3
  egress on small result files. Because the bucket is not Requester Pays, these costs
  fall on the bucket owner (you) — which is correct here, since you own both the bucket
  and the compute.
