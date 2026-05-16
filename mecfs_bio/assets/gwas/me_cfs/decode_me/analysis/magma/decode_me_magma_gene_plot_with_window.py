"""
This task plots the gene-level p-values produced by MAGMA analysis of DecodeME,

A window of 35/10 kb is used to associated variants with genes
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.magma.decode_me_gwas_1_build_37_magma_ensembl_gene_analysis_with_window import (
    DECODE_ME_GWAS_1_MAGMA_ENSEMBL_GENE_ANALYSIS_WITH_WINDOW,
)
from mecfs_bio.assets.reference_data.ensembl_biomart.gene_thesaurus import (
    GENE_THESAURUS,
)
from mecfs_bio.build_system.task.gene_manhattan_plot_task import (
    GeneManhattanPlotTask,
    MagmaGeneSource,
)

DECODE_ME_MAGMA_GENE_PLOT_WITH_WINDOW = GeneManhattanPlotTask.create(
    asset_id="decode_me_magma_gene_plot_with_window",
    source=MagmaGeneSource(
        magma_task=DECODE_ME_GWAS_1_MAGMA_ENSEMBL_GENE_ANALYSIS_WITH_WINDOW,
        gene_thesaurus_task=GENE_THESAURUS,
    ),
)
