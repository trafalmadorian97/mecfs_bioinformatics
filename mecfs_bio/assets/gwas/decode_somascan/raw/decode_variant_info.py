"""
deCODE's variant_info_iceland table (per-allele imputed frequencies in Iceland): the raw text as
published, the deduplicated, sorted parquet, and its clean subset, which decides which deCODE
SomaScan rows are two-allele tests and gives their true effect allele frequency.
"""

from mecfs_bio.assets.gwas.decode_somascan.raw.decode_s3 import DECODE_VARIANT_INFO_REF
from mecfs_bio.build_system.task.decode_database.build_decode_variant_info_task import (
    BuildDecodeCleanVariantInfoTask,
    BuildDecodeVariantInfoTask,
    DownloadDecodeVariantInfoTask,
)

DECODE_VARIANT_INFO_RAW = DownloadDecodeVariantInfoTask.create(
    "decode_variant_info_iceland_raw", DECODE_VARIANT_INFO_REF
)

DECODE_VARIANT_INFO = BuildDecodeVariantInfoTask.create(
    "decode_variant_info_iceland", DECODE_VARIANT_INFO_RAW
)

DECODE_CLEAN_VARIANT_INFO = BuildDecodeCleanVariantInfoTask.create(
    "decode_variant_info_iceland_clean", DECODE_VARIANT_INFO
)
