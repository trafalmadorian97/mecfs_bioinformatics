from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.bmi.pulit_et_al_2018.analysis.pulit_bmi_ppp_rg import \
    PULIT_ET_AL_2018_BMI_PPP_RG_CIS_EXCLUDED
from mecfs_bio.assets.gwas.bmi.pulit_et_al_2018.analysis.pulit_bmi_standard_analysis import \
    PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS
from mecfs_bio.figures.figure_task_list import MULTI_TISSUE_GENE_EXPRESSION_REF, MULTI_TISSUE_CHROMATIN_REF, IMMGEN_REF, \
    CORCES_ATAC_REF, CAHOY_CNS_REF, GTEX_BRAIN_REF
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
            # PULIT_ET_AL_2018_BMI_PPP_RG_CIS_EXCLUDED.display_frame_task
            # PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[MULTI_TISSUE_GENE_EXPRESSION_REF].plot_task_unwrap,
            # PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[MULTI_TISSUE_GENE_EXPRESSION_REF].multiple_testing_task_parquet,
            # PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[MULTI_TISSUE_CHROMATIN_REF].plot_task_unwrap,
            PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[MULTI_TISSUE_CHROMATIN_REF].multiple_testing_task_parquet,
            # PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[
            #     IMMGEN_REF  ].multiple_testing_task_parquet,
            # PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[
            #     CORCES_ATAC_REF ].multiple_testing_task_parquet,
            # PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[
            #     CAHOY_CNS_REF  ].multiple_testing_task_parquet,
            # PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[
            #     GTEX_BRAIN_REF  ].multiple_testing_task_parquet,
        ]
    )
    # DEFAULT_RUNNER.run(PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS.sldsc_tasks.get_terminal_tasks())

if __name__ == "__main__":
    go()