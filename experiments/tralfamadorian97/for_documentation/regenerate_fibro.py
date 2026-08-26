from mecfs_bio.assets.gwas.fibromyalgia.Kerrebijn_et_al.analysis.kerrebeijin_ppp_rg import \
    KERREBEIJIN_ET_AL_PPP_RG_CIS_EXCLUDED
from mecfs_bio.assets.gwas.fibromyalgia.Kerrebijn_et_al.analysis.standard_analysis_kerrebijin_fibro import \
    KERREBIJN_ET_AL_FIBRO_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_ppp_rg import \
    JOHNSTON_ET_AL_PPP_RG_CIS_EXCLUDED
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seronegative.analysis.ra_seronegative_standard_analysis import \
    SERONEGATIVE_RA_STANDARD_ANALYSIS
from mecfs_bio.figures.key_scripts.push_figures import push_figures
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    regenerate_figures(
        [
            # SERONEGATIVE_RA_STANDARD_ANALYSIS.tasks.magma_gene_manhattan_plot_unwrap
            # SERONEGATIVE_RA_STANDARD_ANALYSIS.tasks.heritability_markdown_task_unwrap,
            KERREBEIJIN_ET_AL_PPP_RG_CIS_EXCLUDED.display_frame_task
            # KERREBIJN_ET_AL_FIBRO_STANDARD_ANALYSIS.hba_magma_tasks_unwrap.magma_independent_cluster_plot_svg_unwrap,
            # KERREBIJN_ET_AL_FIBRO_STANDARD_ANALYSIS.hba_magma_tasks_unwrap.extracted_plot_task,

        ]
    )

    push_figures()

if __name__ == '__main__':
    go()