"""
Materialize the HapMap3 CSF hit table (p < 5e-8, 250 kb merge) and summarize it.

Usage: pixi r python experiments/claude/csf_coloc/materialize_csf_hit_table.py
"""

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.csf_pqtl.csf_database.hapmap3.hapmap3_csf_hit_table import (
    HAPMAP_3_CSF_HIT_TABLE_5E8,
)
from mecfs_bio.constants.protein_coloc_constants import HIT_GENE_COL, HIT_PROTEIN_ID_COL


def main() -> None:
    DEFAULT_RUNNER.run(
        [HAPMAP_3_CSF_HIT_TABLE_5E8], incremental_save=True, must_rebuild_transitive=[]
    )
    table = pl.read_parquet(
        DEFAULT_RUNNER.meta_to_path(HAPMAP_3_CSF_HIT_TABLE_5E8.meta)
    )
    print(f"segments {table.height}, aptamers {table[HIT_PROTEIN_ID_COL].n_unique()}")
    with pl.Config(tbl_rows=20):
        print(table.filter(pl.col(HIT_GENE_COL) == "APOE"))


if __name__ == "__main__":
    main()
