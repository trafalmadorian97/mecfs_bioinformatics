"""
This file contains a list of figures to be exported for use in documentation
"""

from mecfs_bio.assets.gwas.ebv_dna.nyeo_et_al.analysis.ebv_dna_standard_analysis import (
    EBV_DNA_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc_plot import (
    CT_LDSC_INITIAL_PLOT,
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
]
