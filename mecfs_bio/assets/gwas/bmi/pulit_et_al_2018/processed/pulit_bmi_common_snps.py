"""
Task to restrict the Pulit et al. (2019) BMI GWAS to variants with MAF >= 0.01.
"""

from mecfs_bio.assets.gwas.bmi.pulit_et_al_2018.raw.pulit_bmi_raw import (
    PULIT_ET_AL_2018_BMI_RAW,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.filter_snps_by_frequency import FilterSNPsFrequencyTask

PULIT_ET_AL_2018_BMI_COMMON_SNPS = FilterSNPsFrequencyTask.create(
    raw_gwas_task=PULIT_ET_AL_2018_BMI_RAW,
    allele_freq_col="Freq_Tested_Allele",
    freq_thresh=0.01,
    id=AssetId("pulit_et_al_2018_bmi_common_snps"),
)
