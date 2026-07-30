from mecfs_bio.assets.gwas.migraine.uk_biobank_2025.raw.uk_biobank_2025_migraine_raw import (
    UK_BIOBANK_2025_MIGRAINE_EUR_DATA_RAW,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.filter_snps_by_frequency import FilterSNPsFrequencyTask

UK_BIOBANK_2025_MIGRAINE_EUR_COMMON_SNPS_TASK = FilterSNPsFrequencyTask.create(
    raw_gwas_task=UK_BIOBANK_2025_MIGRAINE_EUR_DATA_RAW,
    allele_freq_col="effect_allele_frequency",
    freq_thresh=0.01,
    id=AssetId("uk_biobank_2025_migraine_common_snps"),
)
