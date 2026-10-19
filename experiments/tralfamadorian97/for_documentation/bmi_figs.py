from mecfs_bio.assets.gwas.bmi.pulit_et_al_2018.analysis.pulit_bmi_ppp_rg import \
    PULIT_ET_AL_2018_BMI_PPP_RG_CIS_EXCLUDED
from mecfs_bio.assets.gwas.bmi.pulit_et_al_2018.analysis.pulit_bmi_standard_analysis import \
    PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    regenerate_figures(
        [
            # PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS.heritability_markdown_task_unwrap,
            # PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS.ldsc_diagnostic_plot_task_unwrap,
            # PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS.magma_tasks.inner.bar_plot_task,
            # PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS.magma_gene_manhattan_plot_unwrap,
            # PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS.hba_magma_tasks_unwrap.magma_independent_cluster_plot_svg_unwrap,
            # PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS.hba_magma_tasks_unwrap.extracted_plot_task,
            PULIT_ET_AL_2018_BMI_PPP_RG_CIS_EXCLUDED.display_frame_task
        ]
    )

if __name__ == "__main__":
    go()