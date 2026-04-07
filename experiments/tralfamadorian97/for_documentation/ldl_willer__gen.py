from mecfs_bio.assets.gwas.ldl.million_veterans.analysis.mv_ldl_heritability_task import MV_LDL_LDSC_RESULTS_MARKDOWN
from mecfs_bio.assets.gwas.ldl.willer_et_al.analysis.willer_ldl_standard_analysis import \
    WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.mixer.decode_me_univariate_mixer import DECODE_ME_UNIVARIATE_MIXER
from mecfs_bio.assets.gwas.multi_trait.polygenic_overlap.bivariate_mixer.mecfs_pain_bivariate_mixer import \
    MECFS_PAIN_BIVARIATE_MIXER
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.mixer.johnston_et_al_univariate_mixer import \
    JOHNSTON_ET_AL_UNIVARIATE_MIXER
from mecfs_bio.assets.gwas.triglycerides.analysis.triglycide_standard_analysis import \
    WILLER_ET_AL_EUR_TG_STANDARD_ANALYSIS
from mecfs_bio.figures.figure_task_list import MULTI_TISSUE_GENE_EXPRESSION_REF, MULTI_TISSUE_CHROMATIN_REF
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

