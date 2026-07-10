"""
POPs gene prioritization for the Liu et al 2023 IBD meta-analysis.

Runs POPs on top of a POPs-specific MAGMA gene analysis (gene locations derived from
POPs' gene_annot_jun10.txt), as pops.py requires the MAGMA gene set to be a subset of
its own annotation. Used to sanity-check our POPs pipeline against the published IBD
PoPS results from the POPs manuscript.
"""

from mecfs_bio.asset_generator.concrete_pops_asset_generator import (
    concrete_pops_assets_generate,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.magma.liu_et_al_2023_eur_37_magma_pops_gene_analysis import (
    LIU_ET_AL_IBD_2023_EUR_37_MAGMA_POPS_GENE_ANALYSIS,
)

LIU_ET_AL_2023_IBD_POPS = concrete_pops_assets_generate(
    base_name="liu_et_al_2023_ibd_eur_37",
    magma_gene_analysis_task=LIU_ET_AL_IBD_2023_EUR_37_MAGMA_POPS_GENE_ANALYSIS,
)
