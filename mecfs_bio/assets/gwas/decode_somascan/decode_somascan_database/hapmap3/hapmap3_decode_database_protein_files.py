"""
Per-protein deCODE SomaScan files (slim BETA/SE/N over the HapMap3 index, seek index, position
map, significant rows) for every assay of each group.
"""

from mecfs_bio.asset_generator.decode_somascan_protein_asset_generator import (
    generate_decode_somascan_protein_tasks,
)
from mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_index import (
    HAPMAP_3_DECODE_RAW_INDEX,
    HAPMAP_3_DECODE_SMP_INDEX,
    HAPMAP_3_INDEX_NAME,
)

HAPMAP_3_DECODE_SMP_DATABASE = generate_decode_somascan_protein_tasks(
    "smp", HAPMAP_3_DECODE_SMP_INDEX, HAPMAP_3_INDEX_NAME
)
HAPMAP_3_DECODE_RAW_DATABASE = generate_decode_somascan_protein_tasks(
    "raw", HAPMAP_3_DECODE_RAW_INDEX, HAPMAP_3_INDEX_NAME
)
