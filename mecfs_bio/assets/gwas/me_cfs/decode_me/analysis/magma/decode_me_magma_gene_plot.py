from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.magma.decode_me_gwas_1_build_37_magma_ensembl_gene_analysis import (
    DECODE_ME_GWAS_1_MAGMA_ENSEMBL_GENE_ANALYSIS,
)
from mecfs_bio.assets.reference_data.ensembl_biomart.gene_thesaurus import (
    GENE_THESAURUS,
)
from mecfs_bio.build_system.task.gene_manhattan_plot_task import (
    GeneManhattanPlotTask,
    MagmaGeneSource,
)

DECODE_ME_MAGMA_GENE_PLOT = GeneManhattanPlotTask.create(
    asset_id="decode_me_magma_gene_plot",
    source=MagmaGeneSource(
        magma_task=DECODE_ME_GWAS_1_MAGMA_ENSEMBL_GENE_ANALYSIS,
        gene_thesaurus_task=GENE_THESAURUS,
    ),
)
