"""
Asset generator for applying MAGMA gene property analysis using the Human Brain Atlas scRNAseq data as a reference.
"""

from pathlib import PurePath

from attrs import frozen

from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.reference_data.human_brain_atlas.raw.cluster_annotation_term_metadata import (
    CLUSTER_ANNOTATION_TERM_METADATA,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.processed.magma_entrez_gene_locations_data_build_37_unzipped import (
    MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_37_EXTRACTED,
)
from mecfs_bio.assets.reference_data.magma_ld_reference.magma_eur_build_37_1k_genomes_ref_extracted import (
    MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
)
from mecfs_bio.assets.reference_data.magma_specificity_matrices.raw.magma_specificity_matrix_from_hbca_rna_duncan import (
    MAGMA_ENTREZ_SPECIFICITY_MATRIX_HBCA_RNA_DUNCAN,
)
from mecfs_bio.assets.reference_data.research_paper_supplementary_material.duncan_et_al_2025.processed.duncan_et_al_2025_st1_label_columns import (
    DUNCAN_ET_AL_2025_ST1_LABEL_COLS,
)
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameWhiteSpaceSepTextFormat,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.convert_dataframe_to_markdown_task import (
    ConvertDataFrameToMarkdownTask,
)
from mecfs_bio.build_system.task.copy_file_from_directory_task import (
    CopyFileFromDirectoryTask,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.magma.magma_annotate_task import MagmaAnnotateTask
from mecfs_bio.build_system.task.magma.magma_forward_stepwise_select_task import (
    RETAINED_CLUSTERS_COLUMN,
    MagmaForwardStepwiseSelectTask,
)
from mecfs_bio.build_system.task.magma.magma_gene_analysis_task import (
    MagmaGeneAnalysisTask,
)
from mecfs_bio.build_system.task.magma.magma_gene_set_analysis_task import (
    GENE_SET_ANALYSIS_OUTPUT_STEM_NAME,
    MagmaGeneSetAnalysisTask,
    ModelParams,
)
from mecfs_bio.build_system.task.magma.magma_plot_brain_atlas_result_with_stepwise_labels import (
    HBAIndepPlotOptions,
    MAGMAPlotBrainAtlasResultWithStepwiseLabels,
)
from mecfs_bio.build_system.task.magma.magma_snp_location_task import MagmaSNPFileTask
from mecfs_bio.build_system.task.magma.magma_subset_specificity_matrix_using_top_labels import (
    MagmaSubsetSpecificityMatrixWithTopLabels,
)
from mecfs_bio.build_system.task.magma.plot_magma_brain_atlas_result import (
    BRAIN_ATLAS_PLOT_NAME,
    KEY_HBA_ANNOTATION_COLUMNS,
    PlotMagmaBrainAtlasResultTask,
    PlotSettings,
)
from mecfs_bio.build_system.task.multiple_testing_table_task import (
    MultipleTestingTableTask,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe
from mecfs_bio.constants.magma_constants import MAGMA_P_COLUMN


@frozen
class HBAMagmaTasks:
    snp_loc_task: MagmaSNPFileTask
    magma_annotation_task: MagmaAnnotateTask
    magma_gene_analysis_task: MagmaGeneAnalysisTask
    magma_hba_gene_set_task: MagmaGeneSetAnalysisTask
    magma_hba_gene_set_result_extract_task: CopyFileFromDirectoryTask
    magma_hba_result_annotated_task: JoinDataFramesTask
    magma_hba_multiple_testing_task: MultipleTestingTableTask
    magma_hba_result_plot_task: PlotMagmaBrainAtlasResultTask
    extracted_plot_task: CopyFileFromDirectoryTask
    magma_hba_filtered_spec_matrix: MagmaSubsetSpecificityMatrixWithTopLabels | None = (
        None
    )
    magma_hba_conditional_analysis: MagmaGeneSetAnalysisTask | None = None
    magma_hba_forward_select: MagmaForwardStepwiseSelectTask | None = None
    magma_independent_cluster_plot: (
        MAGMAPlotBrainAtlasResultWithStepwiseLabels | None
    ) = None
    independent_clusters_markdown_task: ConvertDataFrameToMarkdownTask | None = None
    magma_independent_clusters_csv: JoinDataFramesTask | None = None

    def terminal_tasks(self) -> list[Task]:
        result: list = [self.extracted_plot_task]
        if self.magma_independent_cluster_plot is not None:
            result += [self.magma_independent_cluster_plot]
        if self.independent_clusters_markdown_task is not None:
            result += [self.independent_clusters_markdown_task]
        return result


def generate_human_brain_atlas_magma_tasks(
    base_name: str,
    gwas_parquet_with_rsids_task: Task,
    sample_size: int,
    plot_settings: PlotSettings,
    include_independent_cluster_plot: bool = False,
    pipes: list[DataProcessingPipe] | None = None,
    hba_indep_plot_options: HBAIndepPlotOptions = HBAIndepPlotOptions(),
) -> HBAMagmaTasks:
    magma_binary_task = MAGMA_1_1_BINARY_EXTRACTED
    gene_loc_file_task = MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_37_EXTRACTED

    p_value_task = MagmaSNPFileTask.create_for_magma_snp_p_value_file_compute_if_needed(
        gwas_parquet_with_rsids_task=gwas_parquet_with_rsids_task,
        asset_id=base_name + "_hba_magma_snp_p_values",
        pipes=pipes,
    )
    snp_loc_task = MagmaSNPFileTask.create_for_magma_snp_pos_file(
        gwas_parquet_with_rsids_task=gwas_parquet_with_rsids_task,
        asset_id=base_name + "_hba_magma_snp_locs",
        pipes=pipes,
    )
    annotations_task = MagmaAnnotateTask.create(
        asset_id=base_name + "_hba_magma_annotations",
        magma_binary_task=magma_binary_task,
        snp_loc_file_task=snp_loc_task,
        gene_loc_file_task=gene_loc_file_task,
        window=(35, 10),  # This choice comes from the Duncan paper
    )
    gene_analysis_task = MagmaGeneAnalysisTask.create(
        asset_id=base_name + "_hba_magma_gene_analysis",
        magma_annotation_task=annotations_task,
        magma_p_value_task=p_value_task,
        magma_binary_task=magma_binary_task,
        magma_ld_ref_task=MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
        ld_ref_file_stem="g1000_eur",
        sample_size=sample_size,
    )

    cell_analysis_task = MagmaGeneSetAnalysisTask.create(
        asset_id=base_name + "_hba_gene_covar",
        magma_gene_analysis_task=gene_analysis_task,
        magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
        gene_set_task=MAGMA_ENTREZ_SPECIFICITY_MATRIX_HBCA_RNA_DUNCAN,
        set_or_covar="covar",
        model_params=ModelParams(direction_covar="greater", condition_hide=[]),
    )

    result_extract_task = CopyFileFromDirectoryTask.create_result_table(
        base_name + "_hba_gene_covar_results_extracted",
        source_directory_task=cell_analysis_task,
        path_inside_directory=PurePath(
            str(GENE_SET_ANALYSIS_OUTPUT_STEM_NAME + ".gsa.out")
        ),
        extension=".txt",
        read_spec=DataFrameReadSpec(DataFrameWhiteSpaceSepTextFormat(comment_code="#")),
    )
    result_annotated_task = JoinDataFramesTask.create_from_result_df(
        asset_id=base_name + "_hba_gene_covar_results_labeled",
        result_df_task=result_extract_task,
        reference_df_task=DUNCAN_ET_AL_2025_ST1_LABEL_COLS,
        how="left",
        left_on=["VARIABLE"],
        right_on=["VARIABLE"],
    )
    multiple_testing_task = MultipleTestingTableTask.create_from_result_table_task(
        alpha=0.01,
        source_task=result_annotated_task,
        method="bonferroni",
        asset_id=base_name + "_hba_gene_covar_results_multiple_testing",
        p_value_column="P",
        apply_filter=False,
    )
    plot_task = PlotMagmaBrainAtlasResultTask.create(
        result_table_task=multiple_testing_task,
        cluster_annotation_task=CLUSTER_ANNOTATION_TERM_METADATA,
        asset_id=base_name + "_hba_atlas_magma_plot",
        plot_settings=plot_settings,
    )
    extracted_plot_task = CopyFileFromDirectoryTask.create_result_table(
        asset_id=base_name + "_hba_magma_plot_extracted",
        source_directory_task=plot_task,
        path_inside_directory=PurePath(BRAIN_ATLAS_PLOT_NAME + ".html"),
        extension=".html",
        read_spec=None,
    )
    if include_independent_cluster_plot:
        magma_hba_filtered_spec_matrix = (
            MagmaSubsetSpecificityMatrixWithTopLabels.create(
                asset_id=base_name + "_hba_magma_filtered_spec_matrix",
                specificity_matrix_task=MAGMA_ENTREZ_SPECIFICITY_MATRIX_HBCA_RNA_DUNCAN,
                magma_gene_covar_analysis_task=cell_analysis_task,
                nominal_sig_level=0.01,
            )
        )
        magma_hba_conditional_analysis = MagmaGeneSetAnalysisTask.create(
            asset_id=base_name + "_hba_magma_conditional_analysis",
            magma_gene_analysis_task=gene_analysis_task,
            magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
            gene_set_task=magma_hba_filtered_spec_matrix,
            set_or_covar="covar",
            model_params=ModelParams(
                direction_covar="greater", condition_hide=[], joint_pairs=True
            ),
        )
        magma_hba_forward_select = MagmaForwardStepwiseSelectTask.create(
            asset_id=base_name + "_hba_magma_forward_stepwise",
            magma_marginal_output_task=cell_analysis_task,
            magma_conditional_output_task=magma_hba_conditional_analysis,
            significance_threshold=0.01,
        )
        magma_independent_clusters_labeled_labeled = (
            JoinDataFramesTask.create_from_result_df(
                asset_id=base_name + "_hba_magma_independent_clusters_labeled",
                result_df_task=magma_hba_forward_select,
                reference_df_task=result_annotated_task,
                how="left",
                right_on=["VARIABLE"],
                left_on=[RETAINED_CLUSTERS_COLUMN],
                out_pipe=SelectColPipe(
                    [RETAINED_CLUSTERS_COLUMN, MAGMA_P_COLUMN]
                    + KEY_HBA_ANNOTATION_COLUMNS
                ),
            )
        )
        independent_clusters_markdown = (
            ConvertDataFrameToMarkdownTask.create_from_result_table_task(
                asset_id=base_name + "_hba_magma_independent_clusters_markdown",
                source_task=magma_independent_clusters_labeled_labeled,
            )
        )

        magma_independent_cluster_plot = (
            MAGMAPlotBrainAtlasResultWithStepwiseLabels.create(
                result_table_task=multiple_testing_task,
                asset_id=base_name + "_hba_magma_independent_cluster_plot",
                stepwise_cluster_list_task=magma_hba_forward_select,
                plot_options=hba_indep_plot_options,
            )
        )
    else:
        magma_hba_filtered_spec_matrix = None
        magma_hba_conditional_analysis = None
        magma_hba_forward_select = None
        magma_independent_cluster_plot = None
        independent_clusters_markdown = None
        magma_independent_clusters_labeled_labeled = None

    return HBAMagmaTasks(
        snp_loc_task=snp_loc_task,
        magma_annotation_task=annotations_task,
        magma_gene_analysis_task=gene_analysis_task,
        magma_hba_gene_set_task=cell_analysis_task,
        magma_hba_gene_set_result_extract_task=result_extract_task,
        magma_hba_result_annotated_task=result_annotated_task,
        magma_hba_multiple_testing_task=multiple_testing_task,
        magma_hba_result_plot_task=plot_task,
        magma_hba_filtered_spec_matrix=magma_hba_filtered_spec_matrix,
        magma_hba_conditional_analysis=magma_hba_conditional_analysis,
        magma_hba_forward_select=magma_hba_forward_select,
        magma_independent_cluster_plot=magma_independent_cluster_plot,
        independent_clusters_markdown_task=independent_clusters_markdown,
        magma_independent_clusters_csv=magma_independent_clusters_labeled_labeled,
        extracted_plot_task=extracted_plot_task,
    )
