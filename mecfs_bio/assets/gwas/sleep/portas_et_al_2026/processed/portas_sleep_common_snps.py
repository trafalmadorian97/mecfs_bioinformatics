"""
Tasks to restrict the Portas et al. (2026) device-measured sleep GWAS to variants with MAF >= 0.01, writing parquet.
"""

from mecfs_bio.assets.gwas.sleep.portas_et_al_2026.raw.portas_sleep_raw import (
    PORTAS_NREM_SLEEP_DURATION_RAW,
    PORTAS_REM_SLEEP_DURATION_RAW,
    PORTAS_SLEEP_DURATION_RAW,
    PORTAS_SLEEP_EFFICIENCY_RAW,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.filter_snps_by_frequency import FilterSNPsFrequencyTask

PORTAS_MAF_THRESHOLD = 0.01

PORTAS_SLEEP_DURATION_COMMON_SNPS = FilterSNPsFrequencyTask.create(
    raw_gwas_task=PORTAS_SLEEP_DURATION_RAW,
    allele_freq_col="effect_allele_frequency",
    freq_thresh=PORTAS_MAF_THRESHOLD,
    id=AssetId("portas_et_al_2026_sleep_duration_common_snps"),
)

PORTAS_SLEEP_EFFICIENCY_COMMON_SNPS = FilterSNPsFrequencyTask.create(
    raw_gwas_task=PORTAS_SLEEP_EFFICIENCY_RAW,
    allele_freq_col="effect_allele_frequency",
    freq_thresh=PORTAS_MAF_THRESHOLD,
    id=AssetId("portas_et_al_2026_sleep_efficiency_common_snps"),
)

PORTAS_REM_SLEEP_DURATION_COMMON_SNPS = FilterSNPsFrequencyTask.create(
    raw_gwas_task=PORTAS_REM_SLEEP_DURATION_RAW,
    allele_freq_col="effect_allele_frequency",
    freq_thresh=PORTAS_MAF_THRESHOLD,
    id=AssetId("portas_et_al_2026_rem_sleep_duration_common_snps"),
)

PORTAS_NREM_SLEEP_DURATION_COMMON_SNPS = FilterSNPsFrequencyTask.create(
    raw_gwas_task=PORTAS_NREM_SLEEP_DURATION_RAW,
    allele_freq_col="effect_allele_frequency",
    freq_thresh=PORTAS_MAF_THRESHOLD,
    id=AssetId("portas_et_al_2026_nrem_sleep_duration_common_snps"),
)
