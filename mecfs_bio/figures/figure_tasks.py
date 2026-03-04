from mecfs_bio.assets.gwas.systemic_lupus_erythematosus.bentham_et_al_2015.analysis_results.bentham_2015_standard_analysis import (
    BENTHAM_LUPUS_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.multi_gwas.genetic_correlation.ct_ldsc_plot import (
    CT_LDSC_INITIAL_PLOT,
)
from mecfs_bio.util.type_related.unwrap import unwrap

ALL_FIGURE_TASKS = [
    CT_LDSC_INITIAL_PLOT,
    unwrap(BENTHAM_LUPUS_STANDARD_ANALYSIS.hba_magma_tasks).extracted_plot_task,
    BENTHAM_LUPUS_STANDARD_ANALYSIS.magma_tasks.inner.bar_plot_task,
]
