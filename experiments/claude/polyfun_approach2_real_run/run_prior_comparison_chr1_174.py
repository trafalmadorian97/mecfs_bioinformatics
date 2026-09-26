"""Build the chr1:174 external-vs-internal prior comparison plots only."""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.prior_comparison_decode_me_37_chr1_174_128_548 import (
    POLYFUN_PRIOR_COMPARISON_CHR1_174,
)

if __name__ == "__main__":
    DEFAULT_RUNNER.run(POLYFUN_PRIOR_COMPARISON_CHR1_174.terminal_tasks())
