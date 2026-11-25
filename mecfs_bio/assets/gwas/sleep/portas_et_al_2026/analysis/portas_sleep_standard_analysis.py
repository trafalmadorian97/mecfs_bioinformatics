"""
Standard analysis of the Portas et al. (2026) device-measured sleep GWAS (sleep duration, sleep efficiency,
REM sleep duration, and non-REM sleep duration).

All four traits are quantitative and share the same 80,013 UK Biobank participants.
"""

from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    StandardAnalysisTaskGroup,
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.sleep.portas_et_al_2026.processed.portas_sleep_common_snps import (
    PORTAS_NREM_SLEEP_DURATION_COMMON_SNPS,
    PORTAS_REM_SLEEP_DURATION_COMMON_SNPS,
    PORTAS_SLEEP_DURATION_COMMON_SNPS,
    PORTAS_SLEEP_EFFICIENCY_COMMON_SNPS,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    QuantPhenotype,
)

PORTAS_SAMPLE_SIZE = 80013

PORTAS_GWAS_CATALOG_HARMONISED_FORMAT = GWASLabColumnSpecifiers(
    chrom="chromosome",
    pos="base_pair_location",
    ea="effect_allele",
    nea="other_allele",
    beta="beta",
    se="standard_error",
    eaf="effect_allele_frequency",
    p="p_value",
    rsid="rsid",
)


def _portas_standard_analysis(
    base_name: str, common_snps_task: Task
) -> StandardAnalysisTaskGroup:
    return concrete_standard_analysis_generator_assume_already_has_rsid(
        base_name=base_name,
        raw_gwas_data_task=common_snps_task,
        fmt=PORTAS_GWAS_CATALOG_HARMONISED_FORMAT,
        sample_size=PORTAS_SAMPLE_SIZE,
        phenotype_info_for_ldsc=QuantPhenotype(),
    )


PORTAS_SLEEP_DURATION_STANDARD_ANALYSIS = _portas_standard_analysis(
    base_name="portas_et_al_2026_sleep_duration",
    common_snps_task=PORTAS_SLEEP_DURATION_COMMON_SNPS,
)

PORTAS_SLEEP_EFFICIENCY_STANDARD_ANALYSIS = _portas_standard_analysis(
    base_name="portas_et_al_2026_sleep_efficiency",
    common_snps_task=PORTAS_SLEEP_EFFICIENCY_COMMON_SNPS,
)

PORTAS_REM_SLEEP_DURATION_STANDARD_ANALYSIS = _portas_standard_analysis(
    base_name="portas_et_al_2026_rem_sleep_duration",
    common_snps_task=PORTAS_REM_SLEEP_DURATION_COMMON_SNPS,
)

PORTAS_NREM_SLEEP_DURATION_STANDARD_ANALYSIS = _portas_standard_analysis(
    base_name="portas_et_al_2026_nrem_sleep_duration",
    common_snps_task=PORTAS_NREM_SLEEP_DURATION_COMMON_SNPS,
)
