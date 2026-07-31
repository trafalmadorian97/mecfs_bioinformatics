"""
Slim per-protein beta/se files over the HapMap3 UKB-PPP database for the European
(discovery) cohort (n = 34,557), European ancestry only.

These assets were built before the sample-size column existed, so include_sample_size is
False (the constant per-protein N comes from a separate ranged-read sample-size table) and
the discovery manifest is passed explicitly to keep them bit-identical to the files already
on disk.
"""

from mecfs_bio.asset_generator.ukbb_ppp_slim_protein_asset_generator import (
    EUR_DISCOVERY_PPP_MANIFEST_PATH,
    generate_ppp_slim_protein_tasks,
)
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.hapmap_3_ppp_index import (
    HAPMAP_3_PPP_DATABASE_INDEX,
)

HAPMAP_3_PPP_DATABASE = generate_ppp_slim_protein_tasks(
    index_task=HAPMAP_3_PPP_DATABASE_INDEX,
    index_name="hapmap_3",
    manifest_path=EUR_DISCOVERY_PPP_MANIFEST_PATH,
    include_sample_size=False,
)
