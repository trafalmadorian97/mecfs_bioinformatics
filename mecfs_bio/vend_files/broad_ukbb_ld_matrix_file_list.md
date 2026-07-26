# broad_ukbb_ld_matrix_file_list.txt

A vendored listing of the Broad Alkes Group's public UK Biobank LD reference
bucket. `ukbb_broad_ld_matrix_generator.py` reads it to discover which genomic
intervals have precomputed LD matrices, then downloads the ones it needs over
plain HTTPS. The listing is checked in so the build does not depend on an S3
call at runtime.

## Regenerating

```bash
aws s3 ls s3://broad-alkesgroup-ukbb-ld/UKBB_LD/ --no-sign-request \
  > mecfs_bio/vend_files/broad_ukbb_ld_matrix_file_list.txt
```

`--no-sign-request` is required: the bucket is public, and without it the CLI
tries to sign with credentials you probably do not have.

The AWS CLI is deliberately **not** a project dependency — it was only ever used
for this one-off listing, so it was removed. Install it temporarily with
`pixi add --feature analysis awscli`, or use any `aws` already on your PATH, and
drop it again afterwards.

## Format the parser depends on

`get_broad_ukbb_ld_matrix_file_info()` parses this positionally, so two
properties of the default `aws s3 ls` output are load-bearing:

- **Whitespace-separated `date time size name`**, one object per line, with the
  filename in the fourth field. Do not add `--human-readable`, `--summarize`, or
  `--recursive` (the last would add path prefixes to the names).
- **Objects in lexicographic order**, which is what S3 returns by default. The
  parser slices `iloc[2:-1]` to drop the two `baselineLF_v2.2.UKB*.tar.gz`
  entries that sort to the top and the `readme_ld.txt` that sorts to the bottom,
  keeping only the `chr<N>_<start>_<end>.{gz,npz,npz2}` rows in between. Sorting
  the file differently silently corrupts the interval table.

As of the last regeneration the listing is 5,557 lines, and the command above
reproduces it byte for byte.
