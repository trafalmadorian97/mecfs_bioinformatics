from mecfs_bio.assets.gwas.schizophrenia.pgc2022.processed.analysis.magma.pgc2022_sch_magma_human_brain_atlas_results_multiple_testing import (
    MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR_MULTIPLE_TESTING,
)
from mecfs_bio.build_system.task.magma.plot_magma_brain_atlas_result import (
    PlotMagmaBrainAtlasResultTask,
)

PGC_2022_SCH_MAGMA_HUMAN_BRAIN_ATLAS_PLOT = PlotMagmaBrainAtlasResultTask.create(
    result_table_task=MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR_MULTIPLE_TESTING,
    asset_id="pgc2022_sch_human_brain_atlas_magma_plot",
)
