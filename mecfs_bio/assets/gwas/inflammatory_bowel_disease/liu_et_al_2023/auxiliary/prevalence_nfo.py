from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
)

LIU_ET_AL_IBD_PREVALENCE_INFO = BinaryPhenotypeSampleInfo(
    sample_prevalence=0.4177,  # 25,042 /(25,042 + 34,915) see Liu et al: https://pmc.ncbi.nlm.nih.gov/articles/PMC10290755/#:~:text=including%2025%2C042%20cases%20and%2034%2C915%20controls%20of%20non%2DFinnish%20European%20(NFE)%20ancestries%20from%20the%20International%20IBD%20Genetics%20Consortium
    estimated_population_prevalence=0.002,  # See Figure 3 of Liu, Zhanju, et al. "Genetic architecture of the inflammatory bowel diseases across East Asian and European ancestries." Nature genetics 55.5 (2023): 796-806.
)
