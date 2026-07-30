from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_standard_analysis import (
    JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.reference_data.ensembl_biomart.gene_thesaurus import (
    GENE_THESAURUS,
)
from mecfs_bio.build_system.task.gene_manhattan_plot_task import (
    GeneManhattanPlotTask,
    MagmaGeneSource,
)

JOHNSTON_MAGMA_GENE_PLOT = GeneManhattanPlotTask.create(
    asset_id="johnston_magma_gene_plot",
    source=MagmaGeneSource(
        magma_task=JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.magma_tasks.inner.gene_analysis_task,
        gene_thesaurus_task=GENE_THESAURUS,
        genome_build="19",
    ),
)
