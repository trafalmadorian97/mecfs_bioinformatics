from mecfs_bio.asset_generator.yu_drg_magma_asset_generator import (
    generate_yu_drg_magma_tasks,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.processed.standard_analysis_sc_pgc_2022 import (
    SCH_PGC_2022_STANDARD_ANALYSIS,
)

PGC_2022_DRG_MAGMA = generate_yu_drg_magma_tasks(
    base_name="pgc_2022",
    gwas_parquet_with_rsids_task=SCH_PGC_2022_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
    sample_size=58749,  # from summary statistics file
)
