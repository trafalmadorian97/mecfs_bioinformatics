"""
Gene set table derived from the MSigDB SQLite database.

One row per human gene set, with all columns from gene_set and gene_set_details,
plus gene_symbols (list[str]) and ncbi_ids (list[int]) aggregated from the
gene_set_gene_symbol / gene_symbol join tables.

SQL query from Claude
"""

from pathlib import PurePath

from mecfs_bio.assets.reference_data.gene_set_data.for_magma.from_gsea_msigdb.extracted.full_msigdb_sqlite_extracted import (
    MSIGDB_SQLLITE_EXTRACTED,
)
from mecfs_bio.build_system.task.sqlite_to_parquet_task import SqliteToParquetTask

_MSIGDB_HUMAN_GENE_SETS_QUERY = """\
SELECT
    gs.id,
    gs.standard_name,
    gs.collection_name,
    gs.tags,
    gs.license_code,
    gsd.systematic_name,
    gsd.description_brief,
    gsd.description_full,
    gsd.exact_source,
    gsd.external_details_URL,
    gsd.source_species_code,
    gsd.primary_namespace_id,
    gsd.second_namespace_id,
    gsd.num_namespaces,
    gsd.publication_id,
    gsd.GEO_id,
    gsd.contributor,
    gsd.contrib_organization,
    gsd.added_in_MSigDB_id,
    gsd.changed_in_MSigDB_id,
    gsd.changed_reason,
    list(sym.symbol ORDER BY gsgs.gene_symbol_id) AS gene_symbols,
    list(TRY_CAST(sym.NCBI_id AS INTEGER) ORDER BY gsgs.gene_symbol_id) AS ncbi_ids
FROM _src.gene_set gs
JOIN _src.gene_set_details gsd ON gsd.gene_set_id = gs.id
JOIN _src.gene_set_gene_symbol gsgs ON gsgs.gene_set_id = gs.id
JOIN _src.gene_symbol sym ON sym.id = gsgs.gene_symbol_id
WHERE gsd.source_species_code = 'HS'
GROUP BY
    gs.id, gs.standard_name, gs.collection_name, gs.tags, gs.license_code,
    gsd.systematic_name, gsd.description_brief, gsd.description_full,
    gsd.exact_source, gsd.external_details_URL, gsd.source_species_code,
    gsd.primary_namespace_id, gsd.second_namespace_id, gsd.num_namespaces,
    gsd.publication_id, gsd.GEO_id, gsd.contributor, gsd.contrib_organization,
    gsd.added_in_MSigDB_id, gsd.changed_in_MSigDB_id, gsd.changed_reason\
"""

MSIGDB_GENE_SETS_PARQUET_FROM_SQLLITE = SqliteToParquetTask.create(
    source_task=MSIGDB_SQLLITE_EXTRACTED,
    asset_id="msigdb_human_gene_sets_table_parquet_from_sqllite",
    query=_MSIGDB_HUMAN_GENE_SETS_QUERY,
    override_subfolder=PurePath("processed"),
)
