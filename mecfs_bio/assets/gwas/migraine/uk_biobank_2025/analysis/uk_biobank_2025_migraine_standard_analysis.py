from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.migraine.uk_biobank_2025.processed.uk_biobank_2025_migraine_common_snps import (
    UK_BIOBANK_2025_MIGRAINE_EUR_COMMON_SNPS_TASK,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
)

UK_BIOBANK_2025_EUR_MIGRAINE_STANDARD_ANALYSIS = concrete_standard_analysis_generator_assume_already_has_rsid(
    base_name="uk_biobank_2025_migraine_eur",
    raw_gwas_data_task=UK_BIOBANK_2025_MIGRAINE_EUR_COMMON_SNPS_TASK,
    fmt=GWASLabColumnSpecifiers(
        chrom="chromosome",
        pos="base_pair_location",
        ea="effect_allele",
        nea="other_allele",
        beta="beta",
        se="standard_error",
        eaf="effect_allele_frequency",
        p="p_value",
        rsid="rsid",
    ),
    sample_size=458440,
    phenotype_info_for_ldsc=BinaryPhenotypeSampleInfo(
        sample_prevalence=25393 / 458440,  # from gwas catalog page
        estimated_population_prevalence=0.104,  # UK prevalence from https://pmc.ncbi.nlm.nih.gov/articles/PMC11753071/
    ),
)
