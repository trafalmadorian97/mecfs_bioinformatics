from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
)

DECODE_ME_PREVALENCE_INFO = BinaryPhenotypeSampleInfo(
    sample_prevalence=0.0566,  # 15,579/(259,909+15,579).  See: https://www.medrxiv.org/content/10.1101/2025.08.06.25333109v1.full-text#:~:text=Our%20primary%20GWAS%2C%20GWAS%2D1%2C%20compared%2015%2C579%20DecodeME%20cases%20with%20259%2C909%20UKB%20controls%20(case%3Acontrol%20ratio%20of%201%3A17)%2C%20across%20all%20autosomes
    estimated_population_prevalence=0.006,  # See Samms and Ponting: https://pmc.ncbi.nlm.nih.gov/articles/PMC12120426/#:~:text=our%20estimated%20ME/CFS%20prevalence%20in%20the%20UK%20then%20rises%20from%20330%2C000%20to%20410%2C000%20(0.6%25)
    total_sample_size=259_909 + 15_579,
)
