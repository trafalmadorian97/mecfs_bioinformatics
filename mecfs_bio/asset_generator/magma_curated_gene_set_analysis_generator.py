from pathlib import PurePath

from mecfs_bio.assets.executable.extracted.magma_binary_extracted import MAGMA_1_1_BINARY_EXTRACTED
from mecfs_bio.assets.reference_data.magma_gene_locations.processed.magma_entrez_gene_locations_data_build_37_unzipped import \
    MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_37_EXTRACTED
from mecfs_bio.assets.reference_data.magma_ld_reference.magma_eur_build_37_1k_genomes_ref_extracted import \
    MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED
from mecfs_bio.assets.reference_data.magma_specificity_matrices.processed.curated_potential_mecfs_gene_sets_specificity_matrix import \
    CURATED_POTENTIAL_MECFS_GENE_SETS_SPECIFICITY_MATRIX
from mecfs_bio.assets.reference_data.magma_specificity_matrices.processed.curated_potential_mecfs_gene_sets_specificity_matrix_reduced import \
    CURATED_POTENTIAL_MECFS_GENE_SETS_SPECIFICITY_MATRIX_REDUCED
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec, \
    DataFrameWhiteSpaceSepTextFormat
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.copy_file_from_directory_task import CopyFileFromDirectoryTask
from mecfs_bio.build_system.task.magma.magma_annotate_task import MagmaAnnotateTask
from mecfs_bio.build_system.task.magma.magma_gene_analysis_task import MagmaGeneAnalysisTask
from mecfs_bio.build_system.task.magma.magma_gene_set_analysis_task import MagmaGeneSetAnalysisTask, ModelParams, \
    GENE_SET_ANALYSIS_OUTPUT_STEM_NAME
from mecfs_bio.build_system.task.magma.magma_snp_location_task import MagmaSNPFileTask
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


def generate_curated_gene_set_analysis_magma_tasks(
        base_name: str,
        gwas_parquet_with_rsids_task: Task,
        sample_size: int,
        pipes: list[DataProcessingPipe] | None = None,
):


    magma_binary_task = MAGMA_1_1_BINARY_EXTRACTED
    gene_loc_file_task = MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_37_EXTRACTED

    p_value_task = MagmaSNPFileTask.create_for_magma_snp_p_value_file_compute_if_needed(
        gwas_parquet_with_rsids_task=gwas_parquet_with_rsids_task,
        asset_id=base_name + "_curated_gene_sets_magma_snp_p_values",
        pipes=pipes,
    )
    snp_loc_task = MagmaSNPFileTask.create_for_magma_snp_pos_file(
        gwas_parquet_with_rsids_task=gwas_parquet_with_rsids_task,
        asset_id=base_name + "_curated_gene_sets_magma_snp_locs",
        pipes=pipes,
    )
    annotations_task = MagmaAnnotateTask.create(
        asset_id=base_name + "_curated_gene_sets_magma_annotations",
        magma_binary_task=magma_binary_task,
        snp_loc_file_task=snp_loc_task,
        gene_loc_file_task=gene_loc_file_task,
        window=(35, 10),  # This choice comes from the Duncan paper
    )
    gene_analysis_task = MagmaGeneAnalysisTask.create(
        asset_id=base_name + "_curated_gene_sets_magma_gene_analysis",
        magma_annotation_task=annotations_task,
        magma_p_value_task=p_value_task,
        magma_binary_task=magma_binary_task,
        magma_ld_ref_task=MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
        ld_ref_file_stem="g1000_eur",
        sample_size=sample_size,
    )
    gene_set_analysis_task_full = MagmaGeneSetAnalysisTask.create(
        asset_id=base_name + "_curated_gene_sets_full_gene_covar",
        magma_gene_analysis_task=gene_analysis_task,
        magma_binary_task=magma_binary_task,
        gene_set_task=CURATED_POTENTIAL_MECFS_GENE_SETS_SPECIFICITY_MATRIX,
        set_or_covar="covar",
        model_params=ModelParams(direction_covar="greater", condition_hide=[]),
    )


    gene_set_analysis_task_reduced = MagmaGeneSetAnalysisTask.create(
        asset_id=base_name + "_curated_gene_sets_reduced_gene_covar",
        magma_gene_analysis_task=gene_analysis_task,
        magma_binary_task=magma_binary_task,
        gene_set_task=CURATED_POTENTIAL_MECFS_GENE_SETS_SPECIFICITY_MATRIX_REDUCED,
        set_or_covar="covar",
        model_params=ModelParams(direction_covar="greater", condition_hide=[]),
    )
    full_result_extraction_task = CopyFileFromDirectoryTask.create_result_table(
        base_name + "_curated_gene_sets_full_covar_results_extracted",
        source_directory_task=gene_set_analysis_task_full,
        path_inside_directory=PurePath(
            str(GENE_SET_ANALYSIS_OUTPUT_STEM_NAME + ".gsa.out")
        ),
        extension=".txt",
        read_spec=DataFrameReadSpec(DataFrameWhiteSpaceSepTextFormat(comment_code="#")),
    )