from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
)

PGC_2022_PREVALENCE_INFO = BinaryPhenotypeSampleInfo(
    sample_prevalence=0.408,  # 53386/(53386 +77258)
    estimated_population_prevalence=0.01,  # See: https://pmc.ncbi.nlm.nih.gov/articles/PMC3327879/#:~:text=Schizophrenia%20is%20a%20severe%20mental,variation%20captured%20by%20common%20SNPs
)
