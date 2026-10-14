from mecfs_bio.assets.gwas.migraine.finngen_r13_migraine_triptan.raw.finngen_r13_migraine_triptan_raw import (
    FINNGEN_R13_MIGRAINE_TRIPTAN_DATA_RAW,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.filter_snps_by_frequency import FilterSNPsFrequencyTask

FINNGEN_R13_MIGRAINE_TRIPTAN_COMMON_SNPS_TASK = FilterSNPsFrequencyTask.create(
    raw_gwas_task=FINNGEN_R13_MIGRAINE_TRIPTAN_DATA_RAW,
    allele_freq_col="af_alt",
    freq_thresh=0.01,
    id=AssetId("finngen_r13_migraine_triptan_common_snps"),
)
