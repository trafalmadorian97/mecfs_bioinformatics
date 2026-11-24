"""
The shared deCODE SomaScan variant index of each group in HapMap3 mode, templated on the first
assay of the group's manifest.
"""

from mecfs_bio.asset_generator.decode_somascan_manifest import (
    load_decode_somascan_proteins,
)
from mecfs_bio.assets.gwas.decode_somascan.raw.decode_variant_info import (
    DECODE_VARIANT_INFO,
)
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.hapmap3_ppp_membership_list import (
    HAPMAP_3_MEMBERSHIP_LIST,
)
from mecfs_bio.build_system.task.decode_database.construct_decode_variant_index_task import (
    ConstructDecodeVariantIndexTask,
)

HAPMAP_3_INDEX_NAME = "hapmap_3"

HAPMAP_3_DECODE_SMP_INDEX = ConstructDecodeVariantIndexTask.create(
    asset_id="decode_somascan_smp_hapmap_3_index",
    template=load_decode_somascan_proteins("smp")[0],
    membership_task=HAPMAP_3_MEMBERSHIP_LIST,
    variant_info_task=DECODE_VARIANT_INFO,
)
HAPMAP_3_DECODE_RAW_INDEX = ConstructDecodeVariantIndexTask.create(
    asset_id="decode_somascan_raw_hapmap_3_index",
    template=load_decode_somascan_proteins("raw")[0],
    membership_task=HAPMAP_3_MEMBERSHIP_LIST,
    variant_info_task=DECODE_VARIANT_INFO,
)
