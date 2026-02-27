from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.magma.decode_me_gwas_1_build_37_magma_ensembl_gene_analysis import (
    DECODE_ME_GWAS_1_MAGMA_ENSEMBL_GENE_ANALYSIS,
)
from mecfs_bio.assets.reference_data.rna_seq_data.gtex.gtex_v10_median_tissue_expression_rna_seq_prep_for_magma import (
    GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA,
)

# from mecfs_bio.assets.reference_data.gene_set_data.for_magma.from_fuma.gtex_v8_tissue_type_gene_sets import (
#     GTEx_V8_TISSUE_EXPRESSION_BASED_GENE_SET_DATA,
# )
from mecfs_bio.build_system.task.magma.magma_gene_set_analysis_task import (
    MagmaGeneSetAnalysisTask,
    ModelParams,
)

MAGMA_DECODE_ME_SPECIFIC_TISSUE_GENE_SET_ANALYSIS = MagmaGeneSetAnalysisTask.create(
    asset_id="decode_me_gwas_1_build_37_magma_ensemble_specific_tissue_gene_covar_analysis",
    magma_gene_analysis_task=DECODE_ME_GWAS_1_MAGMA_ENSEMBL_GENE_ANALYSIS,
    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
    gene_set_task=GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA,
    set_or_covar="covar",
    model_params=ModelParams(direction_covar="greater", condition_hide=["Average"]),
)
