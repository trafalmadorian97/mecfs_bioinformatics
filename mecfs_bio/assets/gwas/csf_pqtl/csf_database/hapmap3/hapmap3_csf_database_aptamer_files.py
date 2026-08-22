"""
Slim per-aptamer beta/se/N files over the HapMap3 CSF pQTL database (Western et al.
2024): one aligned parquet for each of the 7,008 aptamers, in the shared HapMap3
variant index's row order.
"""

from mecfs_bio.asset_generator.csf_slim_aptamer_asset_generator import (
    generate_csf_slim_aptamer_tasks,
)
from mecfs_bio.assets.gwas.csf_pqtl.csf_database.hapmap3.hapmap3_csf_index import (
    HAPMAP_3_CSF_DATABASE_INDEX,
)

HAPMAP_3_CSF_DATABASE = generate_csf_slim_aptamer_tasks(
    index_task=HAPMAP_3_CSF_DATABASE_INDEX,
    index_name="hapmap_3",
)
