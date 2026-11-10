"""
Genome-wide significant (p < 5e-8) hit segments of every UKB-PPP protein over the HapMap3
database for the European (discovery) cohort, merging significant variants within 250 kb.
Used to screen which proteins to colocalize with a trait locus.
"""

from collections.abc import Sequence

from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.eur_discovery_hapmap3_ppp_database_protein_files import (
    HAPMAP_3_PPP_DATABASE,
)
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.hapmap_3_ppp_index import (
    HAPMAP_3_PPP_DATABASE_INDEX,
)
from mecfs_bio.build_system.task.ppp_database.build_slim_protein_parquet_task import (
    BuildSlimProteinParquetTask,
)
from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    HitTableProtein,
    ProteinHitTableTask,
)
from mecfs_bio.constants.protein_coloc_constants import ProteinId


def ppp_hit_table_proteins(
    tasks: Sequence[BuildSlimProteinParquetTask],
) -> tuple[HitTableProtein, ...]:
    """Label each PPP slim protein file by its Olink id and gene."""
    return tuple(
        HitTableProtein(
            task=task,
            protein_id=ProteinId(str(task.protein.oid)),
            gene=task.protein.gene,
        )
        for task in tasks
    )


HAPMAP_3_PPP_HIT_TABLE_5E8 = ProteinHitTableTask.create(
    asset_id="ppp_hit_segments_hapmap_3_p5e-8_gap250kb",
    proteins=ppp_hit_table_proteins(HAPMAP_3_PPP_DATABASE.protein_tasks),
    index_task=HAPMAP_3_PPP_DATABASE_INDEX,
    trait="ukbb_ppp",
    project="ppp_hit_segments",
)
