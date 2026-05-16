from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
)

DECODE_ME_INFECTIOUS_ONSET_PHENOTYPE_INFO = BinaryPhenotypeSampleInfo(
    sample_prevalence=9_738
    / (
        9_738 + 259_909
    ),  # source: https://www.medrxiv.org/content/10.1101/2025.08.06.25333109v1
    estimated_population_prevalence=0.006,  # See Samms and Ponting: https://pmc.ncbi.nlm.nih.gov/articles/PMC12120426/#:~:text=our%20estimated%20ME/CFS%20prevalence%20in%20the%20UK%20then%20rises%20from%20330%2C000%20to%20410%2C000%20(0.6%25)
    total_sample_size=9_738 + 259_909,
)
