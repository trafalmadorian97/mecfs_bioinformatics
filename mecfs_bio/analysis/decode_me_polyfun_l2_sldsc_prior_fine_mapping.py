"""
Script to run the polyfun-vs-uniform explainability fine mapping at the DecodeME
chr15:54.9Mb locus with the DecodeME-specific L2-regularized S-LDSC prior (PolyFun
Approach 2) in place of the precomputed PolyFun prior.

Builds the DecodeME snpvar prior genome-wide (streaming the ~30GB baseline-LF
bundle for its annotation LD-score members on first run), the per-half tau
weights tables that explain it, and the locus's 8-run outer group with its
contrast and plot tasks. Compare against POLYFUN_EXPLAIN_CHR15_54, which uses the
precomputed prior on the same inputs.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.decode_me_l2_sldsc_prior import (
    DECODE_ME_L2_SLDSC_SNPVAR_TABLE,
    DECODE_ME_L2_SLDSC_TAU_EVEN_WEIGHTS,
    DECODE_ME_L2_SLDSC_TAU_ODD_WEIGHTS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr15_54_925_638_l2_sldsc_prior import (
    POLYFUN_EXPLAIN_CHR15_54_L2_SLDSC_PRIOR,
)


def run_decode_me_polyfun_l2_sldsc_prior_fine_mapping():
    """
    Function to build the DecodeME L2-regularized S-LDSC prior and run the
    chr15:54.9Mb explainability fine mapping with it.
    """
    DEFAULT_RUNNER.run(
        [
            DECODE_ME_L2_SLDSC_SNPVAR_TABLE,
            DECODE_ME_L2_SLDSC_TAU_ODD_WEIGHTS,
            DECODE_ME_L2_SLDSC_TAU_EVEN_WEIGHTS,
        ]
        + POLYFUN_EXPLAIN_CHR15_54_L2_SLDSC_PRIOR.terminal_tasks(),
    )


if __name__ == "__main__":
    run_decode_me_polyfun_l2_sldsc_prior_fine_mapping()
