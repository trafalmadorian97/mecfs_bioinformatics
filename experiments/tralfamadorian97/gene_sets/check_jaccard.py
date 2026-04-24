from pathlib import Path

from mecfs_bio.constants.curated_gene_set_collections.curated_mecfs_gene_sets import \
    CURATED_POTENTIAL_MECFS_GENE_SETS_ALL
from mecfs_bio.util.gene_set.msigdb_lookup import gene_set_jaccard_index


def go():
    result = gene_set_jaccard_index(
        parquet_path=Path("assets/base_asset_store/reference_data/gene_set_data/msigdb/extracted/msigdb_human_gene_sets_table_parquet_from_sqllite.parquet.zstd"),
        specs=CURATED_POTENTIAL_MECFS_GENE_SETS_ALL
    )
    print(result)
    import pdb; pdb.set_trace()
    print(result)


if __name__ == '__main__':
    go()