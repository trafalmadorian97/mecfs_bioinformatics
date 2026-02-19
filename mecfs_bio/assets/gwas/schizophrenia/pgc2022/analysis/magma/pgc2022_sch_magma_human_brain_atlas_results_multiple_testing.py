from mecfs_bio.assets.gwas.schizophrenia.pgc2022.analysis.magma.pgc2022_sch_magma_human_brain_atlas_results_labeled import (
    MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR_LABELED,
)
from mecfs_bio.build_system.task.multiple_testing_table_task import (
    MultipleTestingTableTask,
)

MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR_MULTIPLE_TESTING = (
    MultipleTestingTableTask.create_from_result_table_task(
        alpha=0.01,
        source_task=MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR_LABELED,
        method="bonferroni",
        asset_id="pgc2022_sch_human_brain_atlas_gene_covar_results_multiple_testing",
        p_value_column="P",
        apply_filter=False,
    )
)
