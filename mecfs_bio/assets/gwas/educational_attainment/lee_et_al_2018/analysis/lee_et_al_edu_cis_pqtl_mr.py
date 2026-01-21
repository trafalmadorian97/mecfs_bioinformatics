"""
Apply MR using cis-pQTL to the EDU GWAS
"""

from mecfs_bio.asset_generator.mr_with_cis_pqtl_asset_generator import (
    QuantOutcomeConfig,
    mr_cis_pqtl_asset_generator,
)
from mecfs_bio.assets.gwas.educational_attainment.lee_et_al_2018.processed_gwas_data.lee_et_al_magma_task_generator import (
    LEE_ET_AL_2018_COMBINED_MAGMA_TASKS,
)

LEE_ET_AL_EDU_CIS_PQTL_MR = mr_cis_pqtl_asset_generator(
    gwas_dataframe_with_rsids=LEE_ET_AL_2018_COMBINED_MAGMA_TASKS.parquet_file_task,
    base_name="lee_et_al_Edu",
    outcome_config=QuantOutcomeConfig(
        sample_size=257841,
    ),
    steiger_config=None,
)
