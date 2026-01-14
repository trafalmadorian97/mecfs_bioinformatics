from mecfs_bio.assets.gwas.asthma.han_et_al_2022.analysis.asthma_standard_analysis import HAN_ASTHMA_STANDARD_ANALYSIS
from mecfs_bio.assets.reference_data.pqtls.processed.sun_et_al_2023_pqtls_combined_extracted import \
    SUN_ET_AL_2023_COMBINED_PQTLS_EXTRACTED
from mecfs_bio.build_system.task.two_sample_mr_task import TwoSampleMRTask

TwoSampleMRTask.create(
    asset_id="two_sample_mr_asthma_han",
    outcome_data_task=HAN_ASTHMA_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
    exposure_data_task=SUN_ET_AL_2023_COMBINED_PQTLS_EXTRACTED,


)