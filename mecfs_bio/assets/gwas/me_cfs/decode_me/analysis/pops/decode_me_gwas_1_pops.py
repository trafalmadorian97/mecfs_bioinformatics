"""
POPs gene prioritization for the DecodeME ME/CFS GWAS.

Runs POPs on top of the DecodeME build-37 MAGMA ensembl gene analysis. The ensembl
(rather than entrez) MAGMA variant is used because POPs matches genes on Ensembl
gene IDs (its gene_annot_jun10.txt is keyed on ENSGID).
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.magma.decode_me_gwas_1_build_37_magma_ensembl_gene_analysis import (
    DECODE_ME_GWAS_1_MAGMA_ENSEMBL_GENE_ANALYSIS,
)
from mecfs_bio.assets.reference_data.pops.features.pops_features_extracted import (
    POPS_FEATURES_EXTRACTED,
)
from mecfs_bio.assets.reference_data.pops.source.pops_source_extracted import (
    POPS_SOURCE_EXTRACTED,
)
from mecfs_bio.build_system.task_generator.pops_task_generator import PopsTaskGenerator

DECODE_ME_GWAS_1_POPS = PopsTaskGenerator.create(
    base_name="decode_me_gwas_1_build_37",
    pops_source_task=POPS_SOURCE_EXTRACTED,
    raw_features_task=POPS_FEATURES_EXTRACTED,
    magma_gene_analysis_task=DECODE_ME_GWAS_1_MAGMA_ENSEMBL_GENE_ANALYSIS,
)
