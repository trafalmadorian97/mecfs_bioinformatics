from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.prior_comparison_decode_me_37_chr1_174_128_548 import \
    POLYFUN_PRIOR_COMPARISON_CHR1_174
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.susie_explain_decode_me_37_chr1_174_128_548_l2_sldsc_prior import \
    POLYFUN_EXPLAIN_CHR1_174_L2_SLDSC_PRIOR
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    regenerate_figures(
        POLYFUN_EXPLAIN_CHR1_174_L2_SLDSC_PRIOR.terminal_tasks()
    )
    # regenerate_figures(
    #
    # POLYFUN_PRIOR_COMPARISON_CHR1_174.terminal_tasks()
    #
    #                    )

if __name__ == "__main__":
    go()
