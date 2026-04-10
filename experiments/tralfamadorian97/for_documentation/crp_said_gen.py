from mecfs_bio.assets.gwas.c_reactive_protein.said_et_al.analysis.said_crp_standard_analysis import \
    SAID_ET_AL_EUR_CRP_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.ldl.million_veterans.analysis.mv_ldl_heritability_task import MV_LDL_LDSC_RESULTS_MARKDOWN
from mecfs_bio.figures.key_scripts.generate_figures import generate_figures


def go():
    generate_figures([
        SAID_ET_AL_EUR_CRP_STANDARD_ANALYSIS.heritability_markdown_task_unwrap,
        SAID_ET_AL_EUR_CRP_STANDARD_ANALYSIS.magma_tasks.inner.bar_plot_task,

    ])


if __name__ == '__main__':
    go()

