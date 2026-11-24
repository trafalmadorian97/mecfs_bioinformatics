"""
Genome-wide significant (p < 5e-8) hit segments of every deCODE SomaScan assay over the HapMap3
index, merging significant variants within 250 kb, per group. Assays are identified by seq id.
The screen sees HapMap3 variants only, so it misses chrX and loci tagged only by rare variants;
each protein's significant_rows.parquet keeps those rows for a later, denser screen.
"""

from mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_database_protein_files import (
    HAPMAP_3_DECODE_RAW_DATABASE,
    HAPMAP_3_DECODE_SMP_DATABASE,
)
from mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_index import (
    HAPMAP_3_DECODE_RAW_INDEX,
    HAPMAP_3_DECODE_SMP_INDEX,
)
from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    ProteinHitTableTask,
)

HAPMAP_3_DECODE_SMP_HIT_TABLE_5E8 = ProteinHitTableTask.create(
    asset_id="decode_somascan_smp_hapmap_3_hit_table_5e8",
    proteins=HAPMAP_3_DECODE_SMP_DATABASE.hit_table_proteins(),
    index_task=HAPMAP_3_DECODE_SMP_INDEX,
    trait="decode_somascan_smp",
    project="decode_hit_segments",
)
HAPMAP_3_DECODE_RAW_HIT_TABLE_5E8 = ProteinHitTableTask.create(
    asset_id="decode_somascan_raw_hapmap_3_hit_table_5e8",
    proteins=HAPMAP_3_DECODE_RAW_DATABASE.hit_table_proteins(),
    index_task=HAPMAP_3_DECODE_RAW_INDEX,
    trait="decode_somascan_raw",
    project="decode_hit_segments",
)
