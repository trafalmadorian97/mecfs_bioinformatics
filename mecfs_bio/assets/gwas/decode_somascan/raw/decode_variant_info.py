"""
deCODE's variant_info_iceland table (per-allele imputed frequencies in Iceland), deduplicated and
sorted; the source of the true effect allele frequency of every deCODE SomaScan row.
"""

from mecfs_bio.assets.gwas.decode_somascan.raw.decode_s3 import DECODE_VARIANT_INFO_REF
from mecfs_bio.build_system.task.decode_database.build_decode_variant_info_task import (
    BuildDecodeVariantInfoTask,
)

DECODE_VARIANT_INFO = BuildDecodeVariantInfoTask.create(
    "decode_variant_info_iceland", DECODE_VARIANT_INFO_REF
)
