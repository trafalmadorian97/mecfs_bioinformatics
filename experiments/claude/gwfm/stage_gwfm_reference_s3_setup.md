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
s3://<BUCKET>/sbayesrc/reference/Imputed13M/v1/ukbEUR_13M_FullLDM.zip     (~192 GiB)
s3://<BUCKET>/sbayesrc/reference/Imputed13M/v1/ref_b37_1588blocks.pos     (~40 KiB)
s3://<BUCKET>/sbayesrc/reference/Imputed13M/v1/annot_baseline2.2_13M.zip  (~531 MiB)
s3://<BUCKET>/sbayesrc/reference/Imputed13M/v1/gene_map_hg38_hg19.txt     (~4.9 MiB)
```

The marker's consumer prefix (what a later `GctbFineMapTask` recursively copies from)
is `s3://<BUCKET>/sbayesrc/reference/Imputed13M/v1/` — the same folder.

## Step 1: Choose a region

Put the bucket in the **same AWS region you will run GWFM compute in**. The 192 GiB
LD file is copied onto every compute instance; keeping bucket and compute co-located
makes that an intra-region transfer (fast, and no cross-region data-transfer charge).
The bulk of the data movement is bucket -> compute, so the region should follow the
reference data, not your location.

`us-east-1` is a reasonable default (widest instance availability, and it avoids the
`LocationConstraint` quirk in Step 2). Whatever you pick, use it consistently for the
bucket and for the SkyPilot launch region.

## Step 2: Create the bucket

Pick a globally-unique, DNS-compliant name (3-63 chars, lowercase, digits, hyphens;
no underscores). Below, `<BUCKET>` and `<REGION>` are placeholders.

For **us-east-1** (must NOT pass a location constraint):

```
pixi r aws s3api create-bucket --bucket <BUCKET> --region us-east-1
```

For **any other region** (must pass the matching location constraint):

```
pixi r aws s3api create-bucket \
  --bucket <BUCKET> \
  --region <REGION> \
  --create-bucket-configuration LocationConstraint=<REGION>
```

Leave Block Public Access at its default (fully on) — this data is private. No
versioning, lifecycle, or bucket policy is required to stage.

## Step 3: Grant the staging identity permission

Staging uses the standard AWS credential chain (environment variables, shared
credentials file, SSO, or an instance profile) — the code passes boto3 no explicit
credentials. Confirm who you are:

```
pixi r aws sts get-caller-identity
```

That identity needs the permissions the code exercises: `head_object` +
`get_object_attributes` (dedup), a multipart `PutObject` (the large streamed
uploads), and `AbortMultipartUpload` (cleanup if a stream fails mid-upload). This
least-privilege policy, scoped to the one bucket, covers them — attach it to the user
or role (replace `<BUCKET>`):

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
      "Resource": "arn:aws:s3:::<BUCKET>"
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
      "Resource": "arn:aws:s3:::<BUCKET>/*"
    }
  ]
}
```

Note for later: the GWFM compute instances (launched by the SkyPilot executor) will
also need **read** access to this bucket to pull the reference. That is granted
separately via those instances' credentials/instance profile, not by this doc.

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

BUCKET = "<BUCKET>"


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
  s3://<BUCKET>/sbayesrc/reference/Imputed13M/v1/
```

You should see `ukbEUR_13M_FullLDM.zip` (~192 GiB), `annot_baseline2.2_13M.zip`
(~531 MiB), `gene_map_hg38_hg19.txt` (~4.9 MiB), and `ref_b37_1588blocks.pos`
(~40 KiB).

Confirm dedup works by running Step 5 again with `must_rebuild_transitive=[task]`:
every file should now log `already staged` and nothing should re-upload.

The marker asset is written into the local asset store (under the `sbayesrc_gwfm /
ld_reference / Imputed13M/v1` reference path) as `gwfm_reference_marker.json`; its
`s3_prefix` should read `s3://<BUCKET>/sbayesrc/reference/Imputed13M/v1/`.

### Idempotency / resuming a failed run

Staging is idempotent at file granularity. If a run dies partway (e.g. the 60 s
source timeout on the big file), just rerun:

- Files that fully uploaded are detected by the size check and skipped.
- A multipart upload that never completed leaves **no** finished object, so `head`
  returns nothing and that file is re-uploaded from scratch. (Incomplete multipart
  parts can linger and accrue storage cost; a lifecycle rule to abort incomplete
  multipart uploads after e.g. 7 days is a reasonable hygiene add, though not
  required to stage.)

## Cost notes

- **Storage:** One Zone-IA is ~$0.01/GiB-month, so the ~193 GiB bundle is roughly
  $2/month. One Zone-IA bills a 128 KiB minimum object size and a 30-day minimum
  duration — negligible here.
- **Upload transfer:** data *into* S3 is free; you are not charged for the
  gctbhub -> S3 ingest beyond any egress your runner's own network provider bills.
- **Later reads:** bucket -> compute is free when they are in the same region — the
  reason Step 1 co-locates them.
