"""Force-rebuild the chr1:174 prior comparison tables after a code change (the
build cache does not track task code).

Usage: rebuild_prior_comparison_tables.py [label ...]
Rebuilds the tables for the given run-config labels (e.g. l10), or all of them
when none are given.
"""

import sys

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.prior_comparison_decode_me_37_chr1_174_128_548 import (
    POLYFUN_PRIOR_COMPARISON_CHR1_174,
)

if __name__ == "__main__":
    groups = POLYFUN_PRIOR_COMPARISON_CHR1_174.groups_by_label
    labels = sys.argv[1:] or list(groups)
    tables = [groups[label].table for label in labels]
    DEFAULT_RUNNER.run(tables, must_rebuild_transitive=tables)
