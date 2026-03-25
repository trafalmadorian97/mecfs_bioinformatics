"""
This file contains a list of figures to be exported for use in documentation
"""

from mecfs_bio.assets.gwas.brainstem.whole_brainstem.xue_et_al.analysis.xue_whole_brainstem_standard_analysis import (
    XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.ebv_dna.nyeo_et_al.analysis.ebv_dna_standard_analysis import (
    EBV_DNA_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.human_herpesvirus_7_dna.kamitaki_et_al_2025.analysis.kamitaki_et_al_2025_standard_analysis import (
    KAMITAKI_ET_AL_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.mixer.decode_me_univariate_mixer import (
    DECODE_ME_UNIVARIATE_MIXER,
)
from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.ct_ldsc_plot import (
    CT_LDSC_INITIAL_PLOT,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.mixer.johnston_et_al_univariate_mixer import (
    JOHNSTON_ET_AL_UNIVARIATE_MIXER,
)
from mecfs_bio.assets.gwas.systemic_lupus_erythematosus.bentham_et_al_2015.analysis_results.bentham_2015_standard_analysis import (
    BENTHAM_LUPUS_STANDARD_ANALYSIS,
)
from mecfs_bio.util.type_related.unwrap import unwrap

MULTI_TISSUE_CHROMATIN_REF = "multi_tissue_chromatin"
MULTI_TISSUE_GENE_EXPRESSION_REF = "multi_tissue_gene_expression"

ALL_FIGURE_TASKS = [
    CT_LDSC_INITIAL_PLOT,
    unwrap(BENTHAM_LUPUS_STANDARD_ANALYSIS.hba_magma_tasks).extracted_plot_task,
    BENTHAM_LUPUS_STANDARD_ANALYSIS.magma_tasks.inner.bar_plot_task,
    BENTHAM_LUPUS_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[
        MULTI_TISSUE_CHROMATIN_REF
    ].plot_task,
    BENTHAM_LUPUS_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[
        MULTI_TISSUE_GENE_EXPRESSION_REF
    ].plot_task,
    unwrap(EBV_DNA_STANDARD_ANALYSIS.tasks.hba_magma_tasks).extracted_plot_task,
    EBV_DNA_STANDARD_ANALYSIS.tasks.magma_tasks.inner.bar_plot_task,
    EBV_DNA_STANDARD_ANALYSIS.tasks.sldsc_tasks.partitioned_tasks[
        MULTI_TISSUE_CHROMATIN_REF
    ].plot_task,
    DECODE_ME_UNIVARIATE_MIXER.power_plot_task,
    DECODE_ME_UNIVARIATE_MIXER.qq_plot_task,
    DECODE_ME_UNIVARIATE_MIXER.result_markdown_table_task,
    JOHNSTON_ET_AL_UNIVARIATE_MIXER.power_plot_task,
    JOHNSTON_ET_AL_UNIVARIATE_MIXER.qq_plot_task,
    JOHNSTON_ET_AL_UNIVARIATE_MIXER.result_markdown_table_task,
    unwrap(KAMITAKI_ET_AL_STANDARD_ANALYSIS.tasks.hba_magma_tasks).extracted_plot_task,
    unwrap(
        KAMITAKI_ET_AL_STANDARD_ANALYSIS.tasks.hba_magma_tasks
    ).independent_clusters_markdown_task,
    KAMITAKI_ET_AL_STANDARD_ANALYSIS.tasks.magma_tasks.inner.bar_plot_task,
    XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS.magma_tasks.inner.bar_plot_task,
    unwrap(XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS.hba_magma_tasks).extracted_plot_task,
    unwrap(
        unwrap(
            XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS.hba_magma_tasks
        ).independent_clusters_markdown_task
    ),
    unwrap(
        unwrap(
            XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS.hba_magma_tasks
        ).magma_independent_cluster_plot
    ),
    XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[
        MULTI_TISSUE_CHROMATIN_REF
    ].plot_task,
]
