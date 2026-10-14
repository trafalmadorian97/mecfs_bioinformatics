"""
Source: https://www.ebi.ac.uk/gwas/studies/GCST90474065
"""

from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
)

UKBB_ANK_SPOND_PHENOTYPE_INFO = BinaryPhenotypeSampleInfo(
    sample_prevalence=2076 / 458440,
    estimated_population_prevalence=0.0025,  # https://www.nature.com/articles/npjgenmed20168,
    total_sample_size=458440,
)
