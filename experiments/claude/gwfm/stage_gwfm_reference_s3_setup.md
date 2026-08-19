# Creating an S3 bucket to stage the GWFM reference bundle

Step-by-step setup for running `StageGwfmReferenceTask` against real S3, so you can
verify staging (upload + per-file dedup + marker) end to end.

## What the task actually does

`StageGwfmReferenceTask.create(bucket=...)` (in
`mecfs_bio/build_system/task/sbayesrc/stage_gwfm_reference_task.py`):

1. For each file in `GWFM_REFERENCE_BUNDLE`, checks whether it is already in the
   bucket at the canonical prefix (a `head_object` + `get_object_attributes` call).
2. If absent or the size differs, streams the file from its gctbhub `source_url`
   straight into S3 (`urlopen(source_url)` piped into boto3 `upload_fileobj`).
3. Writes a small JSON marker asset recording the version, the S3 prefix consumers
   should read from, and the sorted file list.

Two consequences that shape the setup below:

- **The bytes flow gctbhub -> your machine -> S3.** `execute` runs locally, not on
  a remote instance. Nothing is written to local disk (the upload is streamed), but
  the full ~192 GiB LD file passes through your machine's network. See "Where to run
  it" for why an in-region EC2 box is the practical choice.
- **Uploads use the `ONEZONE_IA` storage class and request an S3-managed SHA-256
  checksum.** Both are per-object options the code sets; no special bucket
  configuration is needed to allow them.

### The files and where they land

The canonical prefix is produced by `gwfm_reference_prefix(version)` =
`sbayesrc/reference/{version}/`, with `version = "Imputed13M/v1"`. So the objects
end up at:

```
s3://mecfs-bio-reference-data/sbayesrc/reference/Imputed13M/v1/ukbEUR_13M_FullLDM.zip     (~192 GiB)
s3://mecfs-bio-reference-data/sbayesrc/reference/Imputed13M/v1/ref_b37_1588blocks.pos     (~40 KiB)
s3://mecfs-bio-reference-data/sbayesrc/reference/Imputed13M/v1/annot_baseline2.2_13M.zip  (~531 MiB)
s3://mecfs-bio-reference-data/sbayesrc/reference/Imputed13M/v1/gene_map_hg38_hg19.txt     (~4.9 MiB)
```

The marker's consumer prefix (what a later `GctbFineMapTask` recursively copies from)
is `s3://mecfs-bio-reference-data/sbayesrc/reference/Imputed13M/v1/` — the same folder.

## Step 1: Choose a region

Put the bucket in the **same AWS region you will run GWFM compute in**. The 192 GiB
LD file is copied onto every compute instance; keeping bucket and compute co-located
makes that an intra-region transfer (fast, and no cross-region data-transfer charge).
The bulk of the data movement is bucket -> compute, so the region should follow the
reference data, not your location.

`us-east-1` is a reasonable default (widest instance availability, and it avoids the
`LocationConstraint` quirk in Step 2). Whatever you pick, use it consistently for the
bucket and for the SkyPilot launch region.

## Step 2: Create the bucket and configure access

Pick a globally-unique, DNS-compliant name (3-63 chars, lowercase, digits, hyphens;
no underscores). Below, `mecfs-bio-reference-data` and `<REGION>` are placeholders.

For **us-east-1** (must NOT pass a location constraint):

```
pixi r aws s3api create-bucket --bucket mecfs-bio-reference-data --region us-east-1
```

For **any other region** (must pass the matching location constraint):

```
pixi r aws s3api create-bucket \
  --bucket mecfs-bio-reference-data \
  --region <REGION> \
  --create-bucket-configuration LocationConstraint=<REGION>
```

No versioning or lifecycle rule is required to stage. The two things that *do* need
configuring are read access for collaborators and who pays for their downloads.

### The access model: readable by any AWS account, Requester Pays

The goal is: **any collaborator can pull the reference to their own EC2 instance
using their own AWS credentials, while you (the owner) pay only for storage — not for
their transfer.** That is exactly what S3 **Requester Pays** does: the downloader pays
their own request and data-transfer costs; the owner keeps paying storage.

Two facts about Requester Pays shape the rest of the setup:

- It **rejects anonymous (unsigned) requests.** Every reader must be an authenticated
  AWS principal and must opt in per request with `--request-payer requester`. This is
  fine here — collaborators launch EC2 with their own credentials, so they are never
  anonymous. You give up only credential-less public access, which this workflow does
  not need.
- Because reads are authenticated but you still want *any* AWS account to be able to
  read, the read grant is a bucket policy with `Principal: "*"`. S3 treats a `"*"`
  policy as "public", so Block Public Access must be relaxed enough to allow it (below).

**Enable Requester Pays:**

```
pixi r aws s3api put-bucket-request-payment \
  --bucket mecfs-bio-reference-data \
  --request-payment-configuration Payer=Requester
```

**Relax Block Public Access enough to attach the read policy.** Keep the two ACL
blocks on (you are not using ACLs); turn off only the two that would block a `"*"`
bucket *policy*:

```
pixi r aws s3api put-public-access-block \
  --bucket mecfs-bio-reference-data \
  --public-access-block-configuration \
  BlockPublicAcls=true,IgnorePublicAcls=true,BlockPublicPolicy=false,RestrictPublicBuckets=false
```

**Attach a read-only bucket policy** granting `GetObject` / `ListBucket` to any AWS
principal (replace `mecfs-bio-reference-data`). Writes are *not* granted here, so the bucket stays
writable only by your staging identity (Step 3):

```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Sid": "GwfmPublicReadReference",
      "Effect": "Allow",
      "Principal": "*",
      "Action": ["s3:GetObject", "s3:ListBucket"],
      "Resource": [
        "arn:aws:s3:::mecfs-bio-reference-data",
        "arn:aws:s3:::mecfs-bio-reference-data/*"
      ]
    }
  ]
}
```

```
pixi r aws s3api put-bucket-policy --bucket mecfs-bio-reference-data --policy file://read-policy.json
```

If you would rather not expose read to *every* AWS account, replace `"Principal": "*"`
with `{"AWS": ["arn:aws:iam::<ACCOUNT_ID>:root", ...]}` listing each collaborator's
account. A named-principal policy is **not** considered public, so you can then leave
Block Public Access fully on and skip the `put-public-access-block` step above — at the
cost of editing the policy whenever a collaborator joins.

## Step 3: Grant the staging identity permission

Staging uses the standard AWS credential chain (environment variables, shared
credentials file, SSO, or an instance profile) — the code passes boto3 no explicit
credentials. Confirm who you are:

```
pixi r aws sts get-caller-identity
```

That identity needs the permissions the code exercises: `head_object` +
`get_object_attributes` (dedup), a multipart `PutObject` (the large streamed
uploads), and `AbortMultipartUpload` (cleanup if a stream fails mid-upload). Save this
least-privilege policy, scoped to the one bucket, as `staging-policy.json`:

```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Sid": "GwfmStagingBucketLevel",
      "Effect": "Allow",
      "Action": [
        "s3:ListBucket",
        "s3:GetBucketLocation",
        "s3:ListBucketMultipartUploads"
      ],
      "Resource": "arn:aws:s3:::mecfs-bio-reference-data"
    },
    {
      "Sid": "GwfmStagingObjectLevel",
      "Effect": "Allow",
      "Action": [
        "s3:PutObject",
        "s3:GetObject",
        "s3:GetObjectAttributes",
        "s3:AbortMultipartUpload",
        "s3:ListMultipartUploadParts"
      ],
      "Resource": "arn:aws:s3:::mecfs-bio-reference-data/*"
    }
  ]
}
```

Now attach it as an **inline policy** on the identity `get-caller-identity` reported.
The right command depends on whether that identity is an IAM **user** or a **role**;
read the `Arn` in the previous output to tell which:

- If the `Arn` looks like `arn:aws:iam::<ACCOUNT>:user/<NAME>` — an IAM user — use
  `put-user-policy` with that `<NAME>`:

  ```
  pixi r aws iam put-user-policy \
    --user-name <NAME> \
    --policy-name GwfmStaging \
    --policy-document file://staging-policy.json
  ```

- If the `Arn` looks like `arn:aws:sts::<ACCOUNT>:assumed-role/<ROLE>/<SESSION>` (an
  assumed role — typical for SSO or an EC2 instance profile), attach to the **role**,
  using the `<ROLE>` segment (not the session suffix):

  ```
  pixi r aws iam put-role-policy \
    --role-name <ROLE> \
    --policy-name GwfmStaging \
    --policy-document file://staging-policy.json
  ```

Two notes:

- Attaching a policy is itself an IAM-privileged action; the identity you run these
  commands as needs `iam:PutUserPolicy` / `iam:PutRolePolicy` (an admin identity does).
  If you administer IAM from a *different* identity than the one that will stage, run
  the attach as the admin, naming the staging user/role.
- If your staging identity is the **account root** or already has broad S3 admin
  access, it can write to the bucket as-is and you can skip this attach entirely — the
  policy above only matters when you want the staging identity locked down to
  least privilege.

Note for later: the GWFM compute instances (launched by the SkyPilot executor) pull
the reference using the read grant from Step 2's bucket policy, under their own
credentials/instance profile. Because the bucket is Requester Pays, those reads pass
`--request-payer requester`; the executor's reference-copy step
(`aws s3 cp --recursive` for `s3_inputs` in `build_sky_task`) already emits that flag.

## Step 4: Where to run it

You can run staging from anywhere with the credentials and network reach, but the
practical choice is an **EC2 instance in the same region as the bucket**:

- The ~192 GiB LD file streams gctbhub -> runner -> S3. On a home/office connection
  that is many hours and is bottlenecked by your uplink. An in-region EC2 box has
  far more bandwidth and writes to S3 intra-region.
- One caveat that makes a stable, fast link matter: `upload_from_url` opens the
  source with a 60-second socket timeout. If gctbhub stalls for longer than that
  between reads, the stream aborts mid-file. Reruns are safe (see Step 6), but a
  faster, steadier link makes an abort far less likely on the big file.

If you use an EC2 box, give it an instance profile carrying the Step 3 policy and
run the repo there; otherwise run locally with your own credentials.

## Step 5: Run StageGwfmReferenceTask

Tasks are executed through the build system's runner. Save this as
`experiments/claude/gwfm/stage_reference.py` (set `BUCKET`), then run it with
`pixi r python experiments/claude/gwfm/stage_reference.py`:

```python
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.build_system.task.sbayesrc.stage_gwfm_reference_task import (
    StageGwfmReferenceTask,
)

BUCKET = "mecfs-bio-reference-data"


def main() -> None:
    task = StageGwfmReferenceTask.create(bucket=BUCKET)
    # A first run has no build trace, so the task executes and stages. To force a
    # re-run after the trace exists, pass must_rebuild_transitive=[task]; the
    # per-file S3 dedup still skips anything already uploaded, so this is cheap.
    DEFAULT_RUNNER.run([task])


if __name__ == "__main__":
    main()
```

The task logs one line per bundle file — `uploaded GWFM reference file` (with the
S3-reported `sha256`) or `GWFM reference file already staged`. Expect the run to be
dominated by the LD file transfer.

## Step 6: Verify

List the staged objects and confirm all four are present with sane sizes:

```
pixi r aws s3 ls --recursive --human-readable \
  s3://mecfs-bio-reference-data/sbayesrc/reference/Imputed13M/v1/
```

You should see `ukbEUR_13M_FullLDM.zip` (~192 GiB), `annot_baseline2.2_13M.zip`
(~531 MiB), `gene_map_hg38_hg19.txt` (~4.9 MiB), and `ref_b37_1588blocks.pos`
(~40 KiB).

Confirm dedup works by running Step 5 again with `must_rebuild_transitive=[task]`:
every file should now log `already staged` and nothing should re-upload.

The marker asset is written into the local asset store (under the `sbayesrc_gwfm /
ld_reference / Imputed13M/v1` reference path) as `gwfm_reference_marker.json`; its
`s3_prefix` should read `s3://mecfs-bio-reference-data/sbayesrc/reference/Imputed13M/v1/`.

### Idempotency / resuming a failed run

Staging is idempotent at file granularity. If a run dies partway (e.g. the 60 s
source timeout on the big file), just rerun:

- Files that fully uploaded are detected by the size check and skipped.
- A multipart upload that never completed leaves **no** finished object, so `head`
  returns nothing and that file is re-uploaded from scratch. Its half-uploaded parts
  are handled as described in "Incomplete multipart uploads" below.

### Incomplete multipart uploads

The big LD file goes up as a multipart upload. If one fails partway, its
already-uploaded parts are stored and billed until the upload is aborted, and they
are invisible to `s3 ls` (only completed objects show there).

boto3 handles the common case for you: s3transfer registers an `AbortMultipartUpload`
cleanup immediately after starting each multipart upload (its
`CreateMultipartUploadTask` does this), so any error raised during the transfer — a
source timeout, an S3 rejection — auto-aborts the upload as the exception unwinds. You
do **not** need to abort those by hand, and the staging code deliberately adds no
abort of its own (it would only duplicate this).

The one gap is a hard stop where no client code runs at all: the process is
`kill -9`'d or OOM-killed, or the box loses power mid-upload. That orphans the parts,
because the abort cleanup never gets to run.

**Check for orphans** (lists in-progress uploads; empty output means none):

```
pixi r aws s3api list-multipart-uploads --bucket mecfs-bio-reference-data
```

**Abort one manually**, using the `Key` and `UploadId` from that listing:

```
pixi r aws s3api abort-multipart-upload \
  --bucket mecfs-bio-reference-data \
  --key sbayesrc/reference/Imputed13M/v1/ukbEUR_13M_FullLDM.zip \
  --upload-id <UploadId>
```

**Auto-clean going forward** with a lifecycle rule that aborts any incomplete
multipart upload older than 7 days. This is the durable safety net: S3 enforces it
server-side regardless of what the client does, so it also covers the hard-kill case
above. Save as `abort-incomplete-mpu-lifecycle.json`:

```json
{
  "Rules": [
    {
      "ID": "AbortIncompleteMultipartUploads",
      "Status": "Enabled",
      "Filter": {"Prefix": ""},
      "AbortIncompleteMultipartUpload": {"DaysAfterInitiation": 7}
    }
  ]
}
```

```
pixi r aws s3api put-bucket-lifecycle-configuration \
  --bucket mecfs-bio-reference-data \
  --lifecycle-configuration file://abort-incomplete-mpu-lifecycle.json
```

Permissions: the Step 3 staging policy already allows the list and abort commands
(`s3:ListBucketMultipartUploads`, `s3:AbortMultipartUpload`,
`s3:ListMultipartUploadParts`). Putting the lifecycle rule needs
`s3:PutLifecycleConfiguration`, which is **not** in that policy — it is a one-time
owner/admin action, so run it as the bucket owner rather than the locked-down staging
identity.

## Cost notes

- **Storage:** One Zone-IA is ~$0.01/GiB-month, so the ~193 GiB bundle is roughly
  $2/month. One Zone-IA bills a 128 KiB minimum object size and a 30-day minimum
  duration — negligible here.
- **Upload transfer:** data *into* S3 is free; you are not charged for the
  gctbhub -> S3 ingest beyond any egress your runner's own network provider bills.
- **Later reads (collaborators pulling the reference):** with Requester Pays enabled
  (Step 2), the *downloader* pays their own request + transfer costs, not you. On top
  of that, bucket -> compute is free when both are in the same region — the reason
  Step 1 co-locates them, so an in-region collaborator pays essentially nothing and you
  pay nothing beyond storage.

## Reading under Requester Pays (already handled in code)

Enabling Requester Pays means every read must pass `--request-payer requester`. The
SkyPilot executor already does this: `build_sky_task` in
`mecfs_bio/build_system/wf/remote_executor/skypilot_remote_executor.py` emits
`aws s3 cp --recursive --request-payer requester <s3_uri> <dest>` for each `s3_inputs`
entry, so GWFM compute reads the reference bucket without a 403. Staging and the Step 6
verification are unaffected — they run under the bucket owner's own identity, and the
owner is never charged as (nor required to declare itself) a requester.
