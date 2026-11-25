"""
Genome-wide significant (p < 5e-8) hit segments of every Western et al. 2024 CSF aptamer over
the HapMap3 CSF database, merging significant variants within 250 kb. Used to screen which
aptamers to colocalize with a trait locus. Aptamers are identified by SomaScan seq id, so the
several aptamers of one gene stay distinct.
"""

from collections.abc import Sequence

from mecfs_bio.assets.gwas.csf_pqtl.csf_database.hapmap3.hapmap3_csf_database_aptamer_files import (
    HAPMAP_3_CSF_DATABASE,
)
from mecfs_bio.assets.gwas.csf_pqtl.csf_database.hapmap3.hapmap3_csf_index import (
    HAPMAP_3_CSF_DATABASE_INDEX,
)
from mecfs_bio.build_system.task.csf_database.build_slim_aptamer_parquet_task import (
    WESTERN_CSF_TRAIT_NAME,
    BuildSlimCsfAptamerParquetTask,
)
from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    HitTableProtein,
    ProteinHitTableTask,
)
from mecfs_bio.constants.protein_coloc_constants import ProteinId


def csf_hit_table_proteins(
    tasks: Sequence[BuildSlimCsfAptamerParquetTask],
) -> tuple[HitTableProtein, ...]:
    """Label each CSF slim aptamer file by its seq id and gene symbol."""
    return tuple(
        HitTableProtein(
            task=task,
            protein_id=ProteinId(task.aptamer.seq_id),
            gene=task.aptamer.entrez_gene_symbol,
        )
        for task in tasks
    )


HAPMAP_3_CSF_HIT_TABLE_5E8 = ProteinHitTableTask.create(
    asset_id="hapmap3_csf_hit_table_5e8",
    proteins=csf_hit_table_proteins(HAPMAP_3_CSF_DATABASE.aptamer_tasks),
    index_task=HAPMAP_3_CSF_DATABASE_INDEX,
    trait=WESTERN_CSF_TRAIT_NAME,
    project="csf_hit_segments",
)
