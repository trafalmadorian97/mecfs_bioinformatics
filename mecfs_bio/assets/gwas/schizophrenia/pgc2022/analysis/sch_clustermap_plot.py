from IPython.core.pylabtools import figsize

from mecfs_bio.assets.gwas.schizophrenia.pgc2022.processed.standard_analysis_sc_pgc_2022 import \
    SCH_PGC_2022_STANDARD_ANALYSIS
from mecfs_bio.assets.reference_data.ensembl_biomart.gene_thesaurus import GENE_THESAURUS
from mecfs_bio.assets.reference_data.rna_seq_data.gtex_v10_median_tissue_expression_rna_seq_prep_for_magma import \
    GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA
from mecfs_bio.build_system.task.magma.magma_clustermap_task import ExpressionMatrixClusterMapTask, TopGenesMode, \
    AllClustersMode, GeneWiseSpecificityNormalization, ClusterWiseZScoreNormalization, VlimForGeneWiseSpec, \
    ExtraSeabornPlotOptions, GeneLabelingSpec, AuxGeneNamingDFSpec
from mecfs_bio.build_system.task.pipes.drop_col_pipe import DropColPipe

SCH_PGC_2022_MAGMA_CLUSTER_MAP=ExpressionMatrixClusterMapTask.create(
    spec_matrix_task=GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA,
    asset_id="pgc2022_sch_magma_clustermap",
    gene_select_mode=TopGenesMode(
        num_genes=300,
        source_task=SCH_PGC_2022_STANDARD_ANALYSIS.magma_tasks.inner.gene_analysis_task,
        p_value_thresh=0.01
    ),
    cluster_select_mode=AllClustersMode(),
    gene_col_name="Gene",
    normalizations=[
        GeneWiseSpecificityNormalization(),
        # ClusterWiseZScoreNormalization()
    ],
    plot_params={
        # "vmin":0,
        # "vmax":0.1
             "cmap" : "vlag",
        "figsize":(20,11),
        "metric":'correlation',
            "dendrogram_ratio":(0.1,0.1)

    },
    vlimspec=VlimForGeneWiseSpec(max_multiple=5),
    # plot_params={
    #     "center" : 0,
    #      "cmap" : "vlag",
    #     "vmin":-2,
    #     "vmax":2,
    # }
extra_plot_options=ExtraSeabornPlotOptions(
    hide_dendogram=True,
    hide_gene_labels=True,
    rotate_cluster_labels=True

),
    gene_labeling_spec=GeneLabelingSpec("Gene stable ID","Gene name"),
    aug_gene_naming_df_spec=AuxGeneNamingDFSpec(GENE_THESAURUS),
    spec_matrix_pipe=DropColPipe(
        # []
        ["Average"]
    )

)