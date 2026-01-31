"""
Similar to  MAGMA_DECODE_ME_SPECIFIC_TISSUE_GENE_SET_ANALYSIS, but averages across brain tissues.

The idea is that gene expression is very correlated among brain tissues, in comparison to correlation of a given brain tissue with a non-brain tissue (see for example, Fig 3A from Urbut 2019). Thus, even when controlling for average expression across all body tissues, multiple non-causal brain tissues may still be significant simply due to their high correlation to a different causal brain tissue.

By controlling for average expression across brain tissues, we can see if expression in any brain tissues is associated with a phenotype, and which is not simply due to general brain-wide expression. This may potentially narrow down which specific tissue may be causal.
"""

from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.magma.decode_me_gwas_1_build_37_magma_ensembl_gene_analysis import (
    DECODE_ME_GWAS_1_MAGMA_ENSEMBL_GENE_ANALYSIS,
)
from mecfs_bio.assets.reference_data.rna_seq_data.gtex_v10_median_tissue_expression_rna_seq_prep_for_magma_brain_average import (
    GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA_BRAIN_AVERAGE,
)
from mecfs_bio.build_system.task.magma.magma_gene_set_analysis_task import (
    MagmaGeneSetAnalysisTask,
    ModelParams,
)

MAGMA_DECODE_ME_SPECIFIC_TISSUE_GENE_SET_ANALYSIS_BRAIN_AVERAGE = MagmaGeneSetAnalysisTask.create(
    asset_id="decode_me_gwas_1_build_37_magma_ensemble_specific_tissue_gene_covar_analysis_brain_average",
    magma_gene_analysis_task=DECODE_ME_GWAS_1_MAGMA_ENSEMBL_GENE_ANALYSIS,
    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
    gene_set_task=GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA_BRAIN_AVERAGE,
    set_or_covar="covar",
    model_params=ModelParams(direction_covar="greater", condition_hide=["Average"]),
)
