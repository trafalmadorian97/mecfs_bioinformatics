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

