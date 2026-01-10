from mecfs_bio.assets.gwas.schizophrenia.pgc2022.analysis.magma.pgc2022_sch_magma_human_brain_atlas import (
    MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.analysis.magma.pgc2022_sch_magma_human_brain_atlas_conditional_analysis import (
    MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_CONDITIONAL_ANALYSIS,
)
from mecfs_bio.build_system.task.magma.magma_forward_stepwise_select_task import (
    MagmaForwardStepwiseSelectTask,
)

MAGMA_PGC2022_HBA_FORWARD_STEPWISE_CLUSTER_SELECT = MagmaForwardStepwiseSelectTask.create(
    asset_id="pgc2022_sch_forward_stepwise_cluster_select",
    magma_marginal_output_task=MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR,
    magma_conditional_output_task=MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_CONDITIONAL_ANALYSIS,
    significance_threshold=0.01,
)
