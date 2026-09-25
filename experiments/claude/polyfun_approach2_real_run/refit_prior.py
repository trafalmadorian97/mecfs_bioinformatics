"""Force a refit of the DecodeME Approach-2 prior and rerun the chr1 and chr15 loci.

The build cache does not track task fields, so a changed alpha grid needs
must_rebuild_transitive to take effect.

Run: pixi r python experiments/claude/polyfun_approach2_real_run/refit_prior.py \
    2>&1 | tee experiments/claude/polyfun_approach2_real_run/run_refined_grid.log
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.decode_me_l2_sldsc_prior import (
    DECODE_ME_L2_SLDSC_SNPVAR,
    DECODE_ME_L2_SLDSC_SNPVAR_TABLE,
    DECODE_ME_L2_SLDSC_TAU_EVEN_WEIGHTS,
    DECODE_ME_L2_SLDSC_TAU_ODD_WEIGHTS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.susie_explain_decode_me_37_chr1_174_128_548_l2_sldsc_prior import (
    POLYFUN_EXPLAIN_CHR1_174_L2_SLDSC_PRIOR,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.susie_explain_decode_me_37_chr15_54_925_638_l2_sldsc_prior import (
    POLYFUN_EXPLAIN_CHR15_54_L2_SLDSC_PRIOR,
)

if __name__ == "__main__":
    loci = [POLYFUN_EXPLAIN_CHR1_174_L2_SLDSC_PRIOR, POLYFUN_EXPLAIN_CHR15_54_L2_SLDSC_PRIOR]
    DEFAULT_RUNNER.run(
        [
            DECODE_ME_L2_SLDSC_SNPVAR_TABLE,
            DECODE_ME_L2_SLDSC_TAU_ODD_WEIGHTS,
            DECODE_ME_L2_SLDSC_TAU_EVEN_WEIGHTS,
        ]
        + [task for locus in loci for task in locus.terminal_tasks()],
        must_rebuild_transitive=[DECODE_ME_L2_SLDSC_SNPVAR],
    )
