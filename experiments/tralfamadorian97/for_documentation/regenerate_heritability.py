from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seronegative.analysis.ra_seronegative_standard_analysis import \
    SERONEGATIVE_RA_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.analysis.ra_seropositive_ppp_rg import \
    SEROPOSITIVE_RA_PPP_RG_TASKS_CIS_EXCLUDED
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.eur_discovery_hapmap_3_heritability_figure_table import \
    HAPMAP_3_PPP_HERITABILITY_FIGURE_TABLE
from mecfs_bio.figures.key_scripts.push_figures import push_figures
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    regenerate_figures(
        [
            # SERONEGATIVE_RA_STANDARD_ANALYSIS.tasks.magma_gene_manhattan_plot_unwrap
            # SERONEGATIVE_RA_STANDARD_ANALYSIS.tasks.heritability_markdown_task_unwrap,
            SEROPOSITIVE_RA_PPP_RG_TASKS_CIS_EXCLUDED.display_frame_task
            # HAPMAP_3_PPP_HERITABILITY_FIGURE_TABLE

        ]
    )
    # push_figures()

if __name__ == '__main__':
    go()