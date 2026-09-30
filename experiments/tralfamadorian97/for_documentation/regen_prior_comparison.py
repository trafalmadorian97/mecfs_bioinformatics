from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.prior_comparison_decode_me_37_chr15_54_925_638 import \
    POLYFUN_PRIOR_COMPARISON_CHR15_54
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.prior_comparison_decode_me_37_chr1_174_128_548 import \
    POLYFUN_PRIOR_COMPARISON_CHR1_174
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.prior_comparison_decode_me_37_chr6_26_215_000 import \
    POLYFUN_PRIOR_COMPARISON_CHR6_26
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.prior_comparison_decode_me_37_chr6_97_505_620 import \
    POLYFUN_PRIOR_COMPARISON_CHR6_97
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.susie_explain_decode_me_37_chr15_54_925_638_l2_sldsc_prior import \
    POLYFUN_EXPLAIN_CHR15_54_L2_SLDSC_PRIOR
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.susie_explain_decode_me_37_chr17_50_237_377_l2_sldsc_prior import \
    POLYFUN_EXPLAIN_CHR17_50_L2_SLDSC_PRIOR
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.susie_explain_decode_me_37_chr1_174_128_548_l2_sldsc_prior import \
    POLYFUN_EXPLAIN_CHR1_174_L2_SLDSC_PRIOR
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.susie_explain_decode_me_37_chr20_47_653_230_l2_sldsc_prior import \
    POLYFUN_EXPLAIN_CHR20_47_L2_SLDSC_PRIOR
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.susie_explain_decode_me_37_chr6_26_215_000_l2_sldsc_prior import \
    POLYFUN_EXPLAIN_CHR6_26_L2_SLDSC_PRIOR
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.susie_explain_decode_me_37_chr6_97_505_620_l2_sldsc_prior import \
    POLYFUN_EXPLAIN_CHR6_97_L2_SLDSC_PRIOR
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    # DEFAULT_RUNNER.run([
    #
    #     # POLYFUN_EXPLAIN_CHR6_26_L2_SLDSC_PRIOR.upset_all_polyfun,
    #     # POLYFUN_EXPLAIN_CHR6_97_L2_SLDSC_PRIOR.upset_all_polyfun,
    #     # POLYFUN_EXPLAIN_CHR17_50_L2_SLDSC_PRIOR.upset_all_polyfun,
    #     # POLYFUN_EXPLAIN_CHR20_47_L2_SLDSC_PRIOR .upset_all_polyfun
    #     POLYFUN_EXPLAIN_CHR15_54_L2_SLDSC_PRIOR.upset_all_polyfun,
    #     POLYFUN_EXPLAIN_CHR15_54_L2_SLDSC_PRIOR.upset_cs50_polyfun,
    #     POLYFUN_PRIOR_COMPARISON_CHR15_54.groups_by_label["l10"].plot_svg,
    #     POLYFUN_PRIOR_COMPARISON_CHR15_54.groups_by_label["l10"].table,
    #
    # ])
    regenerate_figures(
        [
            # POLYFUN_EXPLAIN_CHR15_54_L2_SLDSC_PRIOR.upset_all_polyfun,
            # POLYFUN_EXPLAIN_CHR15_54_L2_SLDSC_PRIOR.upset_cs50_polyfun,
            POLYFUN_PRIOR_COMPARISON_CHR15_54.groups_by_label["l10"].plot_svg,
            # POLYFUN_PRIOR_COMPARISON_CHR15_54.groups_by_label["l10"].table,

        ]
    )
    # regenerate_figures(
    #     [
    #         # POLYFUN_PRIOR_COMPARISON_CHR6_26.
    #         POLYFUN_EXPLAIN_CHR6_26_L2_SLDSC_PRIOR.upset_all_polyfun,
    #         POLYFUN_EXPLAIN_CHR6_97_L2_SLDSC_PRIOR.upset_all_polyfun,
    #         POLYFUN_EXPLAIN_CHR15_54_L2_SLDSC_PRIOR.upset_all_polyfun,
    #
    #         # POLYFUN_EXPLAIN_CHR1_174_L2_SLDSC_PRIOR.upset_all_polyfun,
    #         # POLYFUN_EXPLAIN_CHR1_174_L2_SLDSC_PRIOR.upset_cs50_polyfun,
    #         # POLYFUN_PRIOR_COMPARISON_CHR1_174.groups_by_label["l10"].plot_svg,
    #         # POLYFUN_PRIOR_COMPARISON_CHR1_174.groups_by_label["l10"].table,
    #     ]
    #     # POLYFUN_EXPLAIN_CHR1_174_L2_SLDSC_PRIOR.terminal_tasks()
    # )
    # regenerate_figures(
    #
    # POLYFUN_PRIOR_COMPARISON_CHR1_174.terminal_tasks()
    #
    #                    )

if __name__ == "__main__":
    go()
