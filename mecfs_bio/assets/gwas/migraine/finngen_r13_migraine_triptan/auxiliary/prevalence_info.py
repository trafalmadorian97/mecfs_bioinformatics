from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
)

FINNGEN_R13_MIGRAINE_TRIPTAN_PREVALENCE_INFO = BinaryPhenotypeSampleInfo(
    # Exact case and control counts for endpoint MIGRAINE_TRIPTAN from the FinnGen R13 manifest. See: https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv
    sample_prevalence=56_974 / (443_212 + 56_974),
    # Hautakangas et al. 2022 (Nat Genet), migraine GWAS meta-analysis: "We used a migraine population prevalence of 16%" for liability-scale conversion. See: https://pmc.ncbi.nlm.nih.gov/articles/PMC8837554/
    # Triptan purchase is treated as a proxy for migraine, following Bjornsdottir et al. 2023 (Nat Genet), who added triptan purchasers to migraine cases. See: https://pmc.ncbi.nlm.nih.gov/articles/PMC10632135/
    estimated_population_prevalence=0.16,
    # Exact case and control counts for endpoint MIGRAINE_TRIPTAN from the FinnGen R13 manifest. See: https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv
    total_sample_size=443_212 + 56_974,
)
