from pathlib import Path

from mecfs_bio.constants.curated_gene_set_collections.curated_mecfs_gene_sets import \
    CURATED_POTENTIAL_MECFS_GENE_SETS_ALL
from mecfs_bio.util.gene_set.msigdb_diversity import select_diverse_subset
from rich.pretty import pprint


def calc_diverse_subset_of_gene_sets():
    result = select_diverse_subset(
        gene_sets=CURATED_POTENTIAL_MECFS_GENE_SETS_ALL,
        db_path=Path("assets/base_asset_store/reference_data/gene_set_data/msigdb/extracted/msigdb_human_gene_sets_table_parquet_from_sqllite.parquet.zstd"),
        initial_gene_sets=[],
        target_size=15,
    )
    pprint(result)
    import pdb; pdb.set_trace()
    print(result)


if __name__ == "__main__":
    calc_diverse_subset_of_gene_sets()