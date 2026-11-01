"""
Genome-wide significant (p < 5e-8) hit segments of every UKB-PPP protein over the HapMap3
database for the European (discovery) cohort, merging significant variants within 250 kb.
Used to screen which proteins to colocalize with a trait locus.
"""

from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.eur_discovery_hapmap3_ppp_database_protein_files import (
    HAPMAP_3_PPP_DATABASE,
)
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.hapmap_3_ppp_index import (
    HAPMAP_3_PPP_DATABASE_INDEX,
)
from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    PppProteinHitTableTask,
)

HAPMAP_3_PPP_HIT_TABLE_5E8 = PppProteinHitTableTask.create(
    asset_id="ppp_hit_segments_hapmap_3_p5e-8_gap250kb",
    protein_tasks=HAPMAP_3_PPP_DATABASE.protein_tasks,
    index_task=HAPMAP_3_PPP_DATABASE_INDEX,
)
