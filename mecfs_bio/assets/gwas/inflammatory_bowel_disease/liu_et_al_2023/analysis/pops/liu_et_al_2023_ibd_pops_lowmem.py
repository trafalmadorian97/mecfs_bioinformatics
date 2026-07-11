"""
Low-peak-memory POPs gene prioritization for the Liu et al 2023 IBD meta-analysis.

Identical to liu_et_al_2023_ibd_pops but uses the streaming kernel-ridge
reimplementation (low_mem=True), whose peak memory is independent of the
selected-feature count. IBD selects ~14,564 features uncapped, which OOMs the box
under stock POPs; the low-memory path runs it uncapped and model-identically, so
this is the faithful (uncapped) IBD PoPS result. See liu_et_al_2023_ibd_pops for
the stock (capped) comparison run.
"""

from mecfs_bio.asset_generator.concrete_pops_asset_generator import (
    concrete_pops_assets_generate,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.magma.liu_et_al_2023_eur_37_magma_pops_gene_analysis import (
    LIU_ET_AL_IBD_2023_EUR_37_MAGMA_POPS_GENE_ANALYSIS,
)

LIU_ET_AL_2023_IBD_POPS_LOWMEM = concrete_pops_assets_generate(
    base_name="liu_et_al_2023_ibd_eur_37_lowmem",
    magma_gene_analysis_task=LIU_ET_AL_IBD_2023_EUR_37_MAGMA_POPS_GENE_ANALYSIS,
    low_mem=True,
)
