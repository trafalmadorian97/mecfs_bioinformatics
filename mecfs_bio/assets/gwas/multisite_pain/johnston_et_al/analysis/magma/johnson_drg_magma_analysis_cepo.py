from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_standard_analysis import (
    JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.reference_data.rna_seq_data.yu_drg.processed.yu_drg_cepo_specificity_matrix import (
    YU_DRG_CEPO_SPECIFICITY_MATRIX,
)
from mecfs_bio.build_system.task.magma.magma_gene_set_analysis_task import (
    MagmaGeneSetAnalysisTask,
    ModelParams,
)

MAGMA_JOHNSON_DRG_ANALYSIS_CEPO = MagmaGeneSetAnalysisTask.create(
    asset_id="johnson_build_37_magma_ensemble_drg_gene_covar_analysis_cepo",
    magma_gene_analysis_task=JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.magma_tasks.inner.gene_analysis_task,
    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
    gene_set_task=YU_DRG_CEPO_SPECIFICITY_MATRIX,
    set_or_covar="covar",
    model_params=ModelParams(direction_covar="greater", condition_hide=[]),
)
