from mecfs_bio.assets.reference_data.gene_set_data.for_magma.from_gsea_msigdb.processed.full_msigdb_parquet_from_sqlite import \
    MSIGDB_GENE_SETS_PARQUET_FROM_SQLLITE
from mecfs_bio.build_system.task.magma.prepare_gene_sets_for_magma_task import PrepareGeneSetsForMagmaTask
from mecfs_bio.constants.curated_gene_set_collections.curated_mecfs_gene_sets import \
    CURATED_POTENTIAL_MECFS_GENE_SETS_ALL

CURATED_POTENTIAL_MECFS_GENE_SETS_SPECIFICITY_MATRIX = PrepareGeneSetsForMagmaTask.create(
    asset_id="magma_specificity_matrix_for_curated_potential_mecfs_gene_sets",
    gene_sets=CURATED_POTENTIAL_MECFS_GENE_SETS_ALL,
    parquet_db_task=MSIGDB_GENE_SETS_PARQUET_FROM_SQLLITE
)