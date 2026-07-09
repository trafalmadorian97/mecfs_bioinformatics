"""
POPs gene prioritization for the DecodeME ME/CFS GWAS.

Runs POPs on top of the DecodeME build-37 MAGMA ensembl gene analysis. The ensembl
(rather than entrez) MAGMA variant is used because POPs matches genes on Ensembl
gene IDs (its gene_annot_jun10.txt is keyed on ENSGID).
"""

from mecfs_bio.asset_generator.concrete_pops_asset_generator import (
    concrete_pops_assets_generate,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.magma.decode_me_gwas_1_build_37_magma_ensembl_gene_analysis import (
    DECODE_ME_GWAS_1_MAGMA_ENSEMBL_GENE_ANALYSIS,
)

DECODE_ME_GWAS_1_POPS = concrete_pops_assets_generate(
    base_name="decode_me_gwas_1_build_37",
    magma_gene_analysis_task=DECODE_ME_GWAS_1_MAGMA_ENSEMBL_GENE_ANALYSIS,
)
