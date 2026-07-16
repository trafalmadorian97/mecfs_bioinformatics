from mecfs_bio.asset_generator.ukbb_ppp_slim_protein_asset_generator import (
    generate_ppp_slim_protein_tasks,
)
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.hapmap_3_ppp_index import (
    HAPMAP_3_PPP_DATABASE_INDEX,
)

HAPMAP_3_PPP_DATABASE = generate_ppp_slim_protein_tasks(
    index_task=HAPMAP_3_PPP_DATABASE_INDEX,
    index_name="hapmap_3",
)
