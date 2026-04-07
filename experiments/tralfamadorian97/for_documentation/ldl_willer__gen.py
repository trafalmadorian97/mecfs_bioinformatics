from mecfs_bio.assets.gwas.triglycerides.willer_et_al import \
    WILLER_ET_AL_EUR_TG_STANDARD_ANALYSIS
from mecfs_bio.figures.figure_task_list import MULTI_TISSUE_GENE_EXPRESSION_REF
from mecfs_bio.figures.key_scripts.generate_figures import generate_figures


def go():
    generate_figures([
        # WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS.heritability_markdown_task_unwrap,
        # WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS.magma_tasks.inner.bar_plot_task,
        # WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[
        #     MULTI_TISSUE_GENE_EXPRESSION_REF].plot_task_unwrap,
        #
        # WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[
        #     MULTI_TISSUE_CHROMATIN_REF].plot_task_unwrap,

        WILLER_ET_AL_EUR_TG_STANDARD_ANALYSIS.heritability_markdown_task_unwrap,
        WILLER_ET_AL_EUR_TG_STANDARD_ANALYSIS.magma_tasks.inner.bar_plot_task,
        WILLER_ET_AL_EUR_TG_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[
        MULTI_TISSUE_GENE_EXPRESSION_REF].plot_task_unwrap

    ])


if __name__ == '__main__':
    go()

