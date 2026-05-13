from mecfs_bio.asset_generator.magma_curated_gene_set_analysis_generator import (
    curated_gene_set_analysis_magma_tasks_from_gene_analysis,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_hba_magma_analysis import (
    DECODE_ME_HBA_MAGMA_TASKS,
)

DECODE_ME_CURATED_GENE_SET_ANALYSIS = (
    curated_gene_set_analysis_magma_tasks_from_gene_analysis(
        base_name="decode_me_gsa",
        gene_analysis_task=DECODE_ME_HBA_MAGMA_TASKS.magma_gene_analysis_task,
    )
)
