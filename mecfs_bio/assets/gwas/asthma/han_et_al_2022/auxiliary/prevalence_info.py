from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
)

HAN_ET_AL_ASTHMA_PREVALENCE_INFO = BinaryPhenotypeSampleInfo(
    sample_prevalence=0.1638,  # 64,538 /(64,538 + 329,321)
    estimated_population_prevalence=0.117,  # See Johansson et al., 2019:  https://academic.oup.com/hmg/article/28/23/4022/5540983#:~:text=The%20disease%20prevalence%20in%20the%20Caucasian%20participants%20was%2011.7%25
)
