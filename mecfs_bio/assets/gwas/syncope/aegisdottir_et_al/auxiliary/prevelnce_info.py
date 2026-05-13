from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
)

AEGISDOTTIR_SYNCOPE_PREVALENCE_INFO = BinaryPhenotypeSampleInfo(
    sample_prevalence=56_071
    / 946_861,  # see: Genetic variants associated with syncope implicate neural and autonomic processes (Aegisdottir et al.)
    estimated_population_prevalence=0.3,  # This is the prevalence used in the paper
)
