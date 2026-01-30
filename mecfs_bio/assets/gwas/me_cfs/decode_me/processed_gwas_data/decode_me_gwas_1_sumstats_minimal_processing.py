"""
Task to convert the DECODE ME GWAS-1 summary statistics into a GWASLab Sumstats object, while performing minimal processing.
See https://github.com/Cloufield/gwaslab/blob/1f055b632f62ff69e7a9c8a73155771b0dca52d0/src/gwaslab/g_Sumstats.py#L177
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.filtered_snps_gwas_1 import (
    DECODE_ME_FILTER_SNPS_GWAS_1_TASK,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
)

DECODE_ME_GWAS_1_SUMSTATS_MINIMAL_FILTERING = GWASLabCreateSumstatsTask(
    df_source_task=DECODE_ME_FILTER_SNPS_GWAS_1_TASK,
    asset_id=AssetId("decode_me_gwas_1_sumstats_minimal_processing"),
    basic_check=True,
    genome_build="infer",
)
