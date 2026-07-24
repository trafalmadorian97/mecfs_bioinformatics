from mecfs_bio.assets.gwas.blood_pressure.keaton_et_al_diastolic.analysis.keaton_dbp_ppp_rg import \
    KEATON_DBP_PPP_RG_CIS_EXCLUDED
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    regenerate_figures(
        [
            # SERONEGATIVE_RA_STANDARD_ANALYSIS.tasks.magma_gene_manhattan_plot_unwrap
            # SERONEGATIVE_RA_STANDARD_ANALYSIS.tasks.heritability_markdown_task_unwrap,
            KEATON_DBP_PPP_RG_CIS_EXCLUDED.display_frame_task
        ]
    )
    # push_figures()

if __name__ == '__main__':
    go()
