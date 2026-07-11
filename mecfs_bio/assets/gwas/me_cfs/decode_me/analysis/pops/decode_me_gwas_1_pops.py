"""
POPs gene prioritization for the DecodeME ME/CFS GWAS.

Runs POPs on top of a DecodeME build-37 MAGMA gene analysis. A POPs-specific MAGMA
variant is used (gene locations derived from POPs' gene_annot_jun10.txt) rather than
the standard ensembl one, because pops.py requires every MAGMA gene to be present in
its own annotation; the standard Ensembl v102 gene set contains genes POPs does not
know about and crashes it.
"""

from mecfs_bio.asset_generator.concrete_pops_asset_generator import (
    concrete_pops_assets_generate,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.magma.decode_me_gwas_1_build_37_magma_pops_gene_analysis import (
    DECODE_ME_GWAS_1_MAGMA_POPS_GENE_ANALYSIS,
)

DECODE_ME_GWAS_1_POPS = concrete_pops_assets_generate(
    base_name="decode_me_gwas_1_build_37",
    magma_gene_analysis_task=DECODE_ME_GWAS_1_MAGMA_POPS_GENE_ANALYSIS,
)
