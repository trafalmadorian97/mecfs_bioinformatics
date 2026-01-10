from mecfs_bio.assets.gwas.schizophrenia.pgc2022.analysis.magma.pgc2022_sch_magma_human_brain_atlas_results_extract import (
    MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR_EXTRACTED,
)
from mecfs_bio.assets.reference_data.research_paper_supplementary_material.duncan_et_al_2025.processed.duncan_et_al_2025_st1_label_columns import (
    DUNCAN_ET_AL_2025_ST1_LABEL_COLS,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask

MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR_LABELED = (
    JoinDataFramesTask.create_from_result_df(
        asset_id="pgc2022_sch_human_brain_atlas_gene_covar_results_labeled",
        result_df_task=MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR_EXTRACTED,
        reference_df_task=DUNCAN_ET_AL_2025_ST1_LABEL_COLS,
        how="left",
        left_on=["VARIABLE"],
        right_on=["VARIABLE"],
    )
)
