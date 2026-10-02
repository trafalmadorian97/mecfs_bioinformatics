from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
)

FINNGEN_R13_MIGRAINE_NO_AURA_PREVALENCE_INFO = BinaryPhenotypeSampleInfo(
    # Exact case and control counts for endpoint G6_MIGRAINE_NO_AURA from the FinnGen R13 manifest. See: https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv
    sample_prevalence=10_206 / (366_556 + 10_206),
    # Rasmussen and Olesen 1992 (Cephalalgia), clinical interview of a representative Danish general population sample aged 25-64: lifetime prevalence of migraine without aura was 8%. See: https://doi.org/10.1046/j.1468-2982.1992.1204221.x
    estimated_population_prevalence=0.08,
    # Exact case and control counts for endpoint G6_MIGRAINE_NO_AURA from the FinnGen R13 manifest. See: https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv
    total_sample_size=366_556 + 10_206,
)
