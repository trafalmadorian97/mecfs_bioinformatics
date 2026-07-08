"""
Tasks for running POPs (Polygenic Priority Score) to prioritize genes.

POPs (Weeks et al., Nat Genet 2023; https://github.com/FinucaneLab/pops) predicts
which genes underlie a trait by combining MAGMA gene-level association statistics
with a large set of gene features. It is run as two steps: munging a directory of
raw gene features into chunked matrices, then fitting the POPs model against MAGMA
scores to produce per-gene priority scores.

The POPs code is not distributed as a Python package, so it is fetched as a pinned
GitHub source tarball and invoked as a subprocess. See
mecfs_bio.assets.reference_data.pops for the source and feature assets.
"""
