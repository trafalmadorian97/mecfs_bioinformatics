from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
)

MILLION_VETERANS_ANK_SPOND_PHENOTYPE_INFO = BinaryPhenotypeSampleInfo(
    sample_prevalence=1_637
    / 450_630,  # https://www.ebi.ac.uk/gwas/studies/GCST90476232
    total_sample_size=450_630,
    estimated_population_prevalence=0.0025,  # https://www.nature.com/articles/npjgenmed20168,
)
