from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
)

FINNGEN_R13_MIGRAINE_PREVALENCE_INFO = BinaryPhenotypeSampleInfo(
    # Exact case and control counts for endpoint G6_MIGRAINE from the FinnGen R13 manifest. See: https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv
    sample_prevalence=28_504 / (366_556 + 28_504),
    # Hautakangas et al. 2022 (Nat Genet), migraine GWAS meta-analysis: "We used a migraine population prevalence of 16%" for liability-scale conversion. See: https://pmc.ncbi.nlm.nih.gov/articles/PMC8837554/
    estimated_population_prevalence=0.16,
    # Exact case and control counts for endpoint G6_MIGRAINE from the FinnGen R13 manifest. See: https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv
    total_sample_size=366_556 + 28_504,
)
