"""
MAGMA annotation step for the Liu et al 2023 IBD POPs analysis.

Identical to the standard ensembl annotation except that it uses gene locations
derived from POPs' gene_annot_jun10.txt, so the resulting gene set is a subset of
POPs' annotation (required by pops.py). Reuses the same SNP locations as the ensembl
annotation.
"""

from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.magma.liu_et_al_2023_eur_37_magma_snp_locs import (
    LIU_ET_AL_IBD_EUR_BUILD_37_MAGMA_SNP_LOCS,
)
from mecfs_bio.assets.reference_data.pops.magma_gene_loc.pops_magma_gene_loc import (
    POPS_MAGMA_GENE_LOC,
)
from mecfs_bio.build_system.task.magma.magma_annotate_task import MagmaAnnotateTask

LIU_ET_AL_2023_EUR_IBD_37_POPS_ANNOTATIONS = MagmaAnnotateTask.create(
    asset_id="liu_et_al_2023_eur_37_magma_pops_annotations",
    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
    snp_loc_file_task=LIU_ET_AL_IBD_EUR_BUILD_37_MAGMA_SNP_LOCS,
    gene_loc_file_task=POPS_MAGMA_GENE_LOC,
)
