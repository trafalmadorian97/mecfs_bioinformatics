from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.analysis.magma.pgc_2022_magma_hba_spec_matrix_filtered import (
    PGC2022_HBA_MAGMA_SPEC_MATRIX_FILTERED,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.processed.magma.pgc2022_sch_magna_entrez_gene_analysis import (
    PGC2022_SCH_MAGMA_ENTREZ_GENE_ANALYSIS,
)
from mecfs_bio.build_system.task.magma.magma_gene_set_analysis_task import (
    MagmaGeneSetAnalysisTask,
    ModelParams,
)

MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_CONDITIONAL_ANALYSIS = (
    MagmaGeneSetAnalysisTask.create(
        asset_id="pgc2022_sch_human_brain_atlas_gene_conditional_analysis",
        magma_gene_analysis_task=PGC2022_SCH_MAGMA_ENTREZ_GENE_ANALYSIS,
        magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
        gene_set_task=PGC2022_HBA_MAGMA_SPEC_MATRIX_FILTERED,
        set_or_covar="covar",
        model_params=ModelParams(
            direction_covar="greater", condition_hide=[], joint_pairs=True
        ),
    )
)
