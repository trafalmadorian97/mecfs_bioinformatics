from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
)

FINNGEN_ANK_SPOND_PHENOTYPE_INFO = BinaryPhenotypeSampleInfo(
    sample_prevalence=1462
    / (1462 + 164682),  # https://opengwas.io/datasets/finn-b-M13_ANKYLOSPON#
    total_sample_size=1462
    + 164682,  # https://opengwas.io/datasets/finn-b-M13_ANKYLOSPON#
    estimated_population_prevalence=0.0025,  # https://www.nature.com/articles/npjgenmed20168,
)
