from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.analysis.bellenguez_pp_rg import BELLENGUEZ_PPP_RG_CIS_EXCLUDED
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.analysis.pgc2022_ppp_rg import PGC2022_SCH_PPP_RG_CIS_EXCLUDED
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    regenerate_figures(
        [
            # BELLENGUEZ_PPP_RG_CIS_EXCLUDED.display_frame_task
            PGC2022_SCH_PPP_RG_CIS_EXCLUDED.display_frame_task
        ]
    )

if __name__ == '__main__':
    go()