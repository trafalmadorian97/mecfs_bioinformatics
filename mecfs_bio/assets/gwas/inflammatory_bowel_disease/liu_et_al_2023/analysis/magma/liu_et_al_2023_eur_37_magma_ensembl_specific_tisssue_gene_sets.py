from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.magma.liu_et_al_2023_eur_37_magma_ensembl_gene_analysis import (
    LIU_ET_AL_IBD_2023_EUR_37_MAGMA_ENSEMBL_GENE_ANALYSIS,
)
from mecfs_bio.assets.reference_data.rna_seq_data.gtex_v10_median_tissue_expression_rna_seq_prep_for_magma import (
    GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA,
)
from mecfs_bio.build_system.task.magma.magma_gene_set_analysis_task import (
    MagmaGeneSetAnalysisTask,
    ModelParams,
)

LIU_ET_AL_IBD_EUR_37_SPECIFIC_TISSUE_GENE_SET_ANALYSIS = MagmaGeneSetAnalysisTask.create(
    asset_id="liu_et_al_ibd_37_eur_build_37_magma_ensemble_specific_tissue_gene_covar_analysis",
    magma_gene_analysis_task=LIU_ET_AL_IBD_2023_EUR_37_MAGMA_ENSEMBL_GENE_ANALYSIS,
    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
    gene_set_task=GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA,
    set_or_covar="covar",
    model_params=ModelParams(direction_covar="greater", condition_hide=["Average"]),
)
