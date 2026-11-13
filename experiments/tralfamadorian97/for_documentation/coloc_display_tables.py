from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.ldl.willer_et_al.analysis.willer_ldl_csf_coloc import WILLER_LDL_CSF_COLOC
from mecfs_bio.assets.gwas.ldl.willer_et_al.analysis.willer_ldl_ppp_coloc import WILLER_LDL_PPP_COLOC
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_ppp_coloc import DECODE_ME_PPP_COLOC, \
    DECODE_ME_PPP_COLOC_EXPLORATORY, DECODE_ME_PPP_COLOC_EXPLORATORY_2, DECODE_ME_PPP_COLOC_EXPLORATORY_2_STACK_PLOTS, \
    DECODE_ME_PPP_COLOC_EXPLORATORY_2_PROTEIN_INFO
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.analysis.ra_seropositive_ppp_coloc import \
    RA_SEROPOSITIVE_PPP_COLOC
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
   # DEFAULT_RUNNER.run(
   #     [
   #         # DECODE_ME_PPP_COLOC.display_task
   #         DECODE_ME_PPP_COLOC_EXPLORATORY.display_task
   #      ],
   #     must_rebuild_transitive=[
   #         # DECODE_ME_PPP_COLOC.display_task
   #                              ],
   # )
   regenerate_figures(
       [
           # DECODE_ME_PPP_COLOC_EXPLORATORY_2_STACK_PLOTS.png_task,
           DECODE_ME_PPP_COLOC_EXPLORATORY_2_PROTEIN_INFO,
       ]
   )

if __name__ == '__main__':
    go()