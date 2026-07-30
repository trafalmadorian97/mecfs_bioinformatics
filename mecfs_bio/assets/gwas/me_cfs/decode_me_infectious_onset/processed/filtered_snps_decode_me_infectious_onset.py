from mecfs_bio.assets.gwas.me_cfs.decode_me.raw_gwas_data.decode_me_quality_control_snps import (
    DECODE_ME_QC_SNPS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_infectious_onset.raw.decode_me_infectious_onset_raw import (
    DECODE_ME_GWAS_INFECTIOUS_ONSET_TASK,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.filter_snps_task import FilterSNPsTask

DECODE_ME_INFECTIOUS_ONSET_FILTER_SNPS_TASK = FilterSNPsTask.create(
    raw_gwas_task=DECODE_ME_GWAS_INFECTIOUS_ONSET_TASK,
    snp_list_task=DECODE_ME_QC_SNPS,
    id=AssetId("decode_me_infectious_onset_filtered_for_quality_control"),
)
