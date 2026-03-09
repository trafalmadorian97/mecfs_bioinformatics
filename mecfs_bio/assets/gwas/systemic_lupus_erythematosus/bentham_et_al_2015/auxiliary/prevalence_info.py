from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
)

BENTHAM_LUPUS_PREVALENCE_INFO = BinaryPhenotypeSampleInfo(
    sample_prevalence=0.3645,  # https://www.ebi.ac.uk/gwas/studies/GCST003156
    estimated_population_prevalence=0.001,  # rough estimate, based on the wide ranged reported in https://pmc.ncbi.nlm.nih.gov/articles/PMC8982275/pdf/nihms-1780632.pdf
)
