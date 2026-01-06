from attrs import frozen

from mecfs_bio.assets.executable.extracted.magma_binary_extracted import MAGMA_1_1_BINARY_EXTRACTED
from mecfs_bio.assets.reference_data.magma_gene_locations.processed.magma_entrez_gene_locations_data_build_37_unzipped import \
    MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_37_EXTRACTED
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.copy_file_from_directory_task import CopyFileFromDirectoryTask
from mecfs_bio.build_system.task.fdr_multiple_testing_table_task import MultipleTestingTableTask
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.magma.magma_annotate_task import MagmaAnnotateTask
from mecfs_bio.build_system.task.magma.magma_gene_analysis_task import MagmaGeneAnalysisTask
from mecfs_bio.build_system.task.magma.magma_gene_set_analysis_task import MagmaGeneSetAnalysisTask
from mecfs_bio.build_system.task.magma.magma_snp_location_task import MagmaSNPFileTask
from mecfs_bio.build_system.task.magma.plot_magma_brain_atlas_result import PlotMagmaBrainAtlasResultTask


@frozen
class HBAMagmaTasks:
    snp_loc_task:MagmaSNPFileTask
    magma_annotation_task:MagmaAnnotateTask
    magma_gene_analysis_task:MagmaGeneAnalysisTask
    magma_hba_gene_set_task:MagmaGeneSetAnalysisTask
    magma_hba_gene_set_result_extract_task:CopyFileFromDirectoryTask
    magma_hba_result_annotated_task:JoinDataFramesTask
    magma_hba_multiple_testing_task:MultipleTestingTableTask
    magma_hba_result_plot_task:PlotMagmaBrainAtlasResultTask


def generate_human_brain_atlas_magma_tasks(
        base_name:str,
gwas_parquet_with_rsids_task:Task
):

    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED
    snp_loc_task = MagmaSNPFileTask.create_for_magma_snp_pos_file(
        gwas_parquet_with_rsids_task=gwas_parquet_with_rsids_task,
        asset_id=base_name + "hba_magma_snp_locs",
    )
    gene_loc_file_task = MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_37_EXTRACTED

    p_value_task = (
        MagmaSNPFileTask.create_for_magma_snp_p_value_file_compute_if_needed(
            gwas_parquet_with_rsids_task=gwas_parquet_with_rsids_task,
            asset_id=base_name + "_magma_snp_p_values",
        )
    )
    snp_loc_task = MagmaSNPFileTask.create_for_magma_snp_pos_file(
        gwas_parquet_with_rsids_task=gwas_parquet_with_rsids_task,
        asset_id=base_name + "_magma_snp_locs",
    )
    annotations_task = MagmaAnnotateTask.create(
        asset_id=base_name + "_magma_annotations",
        magma_binary_task=magma_binary_task,
        snp_loc_file_task=snp_loc_task,
        gene_loc_file_task=gene_loc_file_task,
        window=(35, 10),
    )
    gene_analysis_task = MagmaGeneAnalysisTask.create(
        asset_id=base_name + "_magma_gene_analysis",
        magma_annotation_task=annotations_task,
        magma_p_value_task=p_value_task,
        magma_binary_task=magma_binary_task,
        magma_ld_ref_task=magma_ld_ref_task,
        ld_ref_file_stem=ld_ref_file_stem,
        sample_size=sample_size,
    )

