from mecfs_bio.assets.gwas.schizophrenia.pgc2022.processed.standard_analysis_sc_pgc_2022 import \
    SCH_PGC_2022_STANDARD_ANALYSIS
from mecfs_bio.assets.reference_data.rna_seq_data.gtex_v10_median_tissue_expression_rna_seq_prep_for_magma import \
    GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA
from mecfs_bio.build_system.task.magma.magma_expression_matrix_task import ExpressionMatrixClusterMapTask, TopGenesMode, \
    AllClustersMode
SCH_PGC_2022_MAGMA_CLUSTER_MAP=ExpressionMatrixClusterMapTask.create(
    spec_matrix_task=GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA,
    asset_id="pgc2022_sch_magma_clustermap",
    gene_select_mode=TopGenesMode(
        num_genes=30,
        source_task=SCH_PGC_2022_STANDARD_ANALYSIS.magma_tasks.inner.gene_analysis_task,
        p_value_thresh=0.01
    ),
    cluster_select_mode=AllClustersMode(),
    gene_col_name="Gene"
)