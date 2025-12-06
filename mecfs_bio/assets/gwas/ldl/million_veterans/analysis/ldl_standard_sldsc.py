from mecfs_bio.asset_generation.sldsc_generator import standard_sldsc_task_generator
from mecfs_bio.assets.gwas.ldl.million_veterans.processed_gwas_data.million_veterans_ldl_eur_magma_task_generator import \
    MILLION_VETERANS_EUR_LDL_MAGMA_TASKS

LDL_STANDARD_SLDSC_TASK_GROUP =standard_sldsc_task_generator(
    sumstats_task=MILLION_VETERANS_EUR_LDL_MAGMA_TASKS.sumstats_task,
    base_name="million_veterans_eur_ldl",
)