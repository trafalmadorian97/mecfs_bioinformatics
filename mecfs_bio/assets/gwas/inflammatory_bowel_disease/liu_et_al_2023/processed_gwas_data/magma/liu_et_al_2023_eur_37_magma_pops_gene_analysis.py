"""
MAGMA gene analysis for the Liu et al 2023 IBD POPs analysis.

Uses the POPs-gene-set annotation (gene locations derived from gene_annot_jun10.txt)
so the output gene set is a subset of POPs' annotation, which pops.py requires. Every
other setting matches the standard ensembl gene analysis.
"""

from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.magma.liu_et_al_2023_eur_37_magma_pops_annotations import (
    LIU_ET_AL_2023_EUR_IBD_37_POPS_ANNOTATIONS,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.magma.liu_et_al_2023_eur_ibd_37_p_values import (
    LIU_ET_AL_2023_IBD_37_MAGMA_SNP_P_VALUES,
)
from mecfs_bio.assets.reference_data.magma_ld_reference.magma_eur_build_37_1k_genomes_ref_extracted import (
    MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
)
from mecfs_bio.build_system.task.magma.magma_gene_analysis_task import (
    MagmaGeneAnalysisTask,
)

LIU_ET_AL_IBD_2023_EUR_37_MAGMA_POPS_GENE_ANALYSIS = MagmaGeneAnalysisTask.create(
    asset_id="liu_et_al_ibd_2023_eur_37_magma_pops_gene_analysis",
    magma_annotation_task=LIU_ET_AL_2023_EUR_IBD_37_POPS_ANNOTATIONS,
    magma_p_value_task=LIU_ET_AL_2023_IBD_37_MAGMA_SNP_P_VALUES,
    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
    magma_ld_ref_task=MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
    ld_ref_file_stem="g1000_eur",
    sample_size=59957,  # from paper. total number of european subjects
)
