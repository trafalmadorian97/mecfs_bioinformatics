"""
Slim per-protein beta/se (plus N) files over the HapMap3 UKB-PPP database for the full
European (combined) cohort (discovery + replication, n = 52,363).

Parallel to eur_discovery_hapmap3_ppp_database_protein_files, but drawn from the Combined
Synapse folder (the eur_combined manifest, the default) and aligned onto the same shared
HapMap3 variant index. The index_name carries the cohort so the assets do not collide with
the discovery ones. include_sample_size defaults on, so N is stored alongside beta/se.
"""

from mecfs_bio.asset_generator.ukbb_ppp_slim_protein_asset_generator import (
    generate_ppp_slim_protein_tasks,
)
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.hapmap_3_ppp_index import (
    HAPMAP_3_PPP_DATABASE_INDEX,
)

HAPMAP_3_PPP_DATABASE_EUR_COMBINED = generate_ppp_slim_protein_tasks(
    index_task=HAPMAP_3_PPP_DATABASE_INDEX,
    index_name="hapmap_3_eur_combined",
)
