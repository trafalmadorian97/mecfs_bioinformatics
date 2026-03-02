from attrs import frozen

from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.raw.magma_ensembl_gene_location_reference_data_build_37 import (
    MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
)
from mecfs_bio.assets.reference_data.magma_ld_reference.magma_eur_build_37_1k_genomes_ref_extracted import (
    MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
)
from mecfs_bio.assets.reference_data.rna_seq_data.yu_drg.processed.yu_drg_cepo_specificity_matrix import (
    YU_DRG_CEPO_SPECIFICITY_MATRIX,
)
from mecfs_bio.assets.reference_data.rna_seq_data.yu_drg.processed.yu_drg_frac_specificity_matrix import (
    YU_DRG_FRAC_SPECIFICITY_MATRIX,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.magma.magma_annotate_task import MagmaAnnotateTask
from mecfs_bio.build_system.task.magma.magma_gene_analysis_task import (
    MagmaGeneAnalysisTask,
)
from mecfs_bio.build_system.task.magma.magma_gene_set_analysis_task import (
    MagmaGeneSetAnalysisTask,
    ModelParams,
)
from mecfs_bio.build_system.task.magma.magma_plot_gene_set_result import (
    MAGMAPlotGeneSetResult,
)
from mecfs_bio.build_system.task.magma.magma_snp_location_task import MagmaSNPFileTask
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class YuDRGMAGMA:
    cepo_plot: Task
    frac_plot: Task

    def terminal_tasks(self) -> list[Task]:
        return [self.cepo_plot, self.frac_plot]


def generate_yu_drg_magma_tasks(
    base_name: str,
    gwas_parquet_with_rsids_task: Task,
    sample_size: int,
    pipes: list[DataProcessingPipe] | None = None,
):
    magma_binary_task = MAGMA_1_1_BINARY_EXTRACTED
    gene_loc_task = MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW

    p_value_task = MagmaSNPFileTask.create_for_magma_snp_p_value_file_compute_if_needed(
        gwas_parquet_with_rsids_task=gwas_parquet_with_rsids_task,
        asset_id=base_name + "_drg_magma_snp_p_values",
        pipes=pipes,
    )
    snp_loc_task = MagmaSNPFileTask.create_for_magma_snp_pos_file(
        gwas_parquet_with_rsids_task=gwas_parquet_with_rsids_task,
        asset_id=base_name + "_drg_magma_snp_locs",
        pipes=pipes,
    )
    annotations_task = MagmaAnnotateTask.create(
        asset_id=base_name + "_drg_magma_annotations",
        magma_binary_task=magma_binary_task,
        snp_loc_file_task=snp_loc_task,
        gene_loc_file_task=gene_loc_task,
        window=(35, 10),  # This choice comes from the Duncan paper
    )
    gene_analysis_task = MagmaGeneAnalysisTask.create(
        asset_id=base_name + "_drg_magma_gene_analysis",
        magma_annotation_task=annotations_task,
        magma_p_value_task=p_value_task,
        magma_binary_task=magma_binary_task,
        magma_ld_ref_task=MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
        ld_ref_file_stem="g1000_eur",
        sample_size=sample_size,
    )

    cell_analysis_task_cepo = MagmaGeneSetAnalysisTask.create(
        asset_id=base_name + "_drg_cell_analysis_cepo",
        magma_gene_analysis_task=gene_analysis_task,
        magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
        gene_set_task=YU_DRG_CEPO_SPECIFICITY_MATRIX,
        set_or_covar="covar",
        model_params=ModelParams(direction_covar="greater", condition_hide=[]),
    )

    cell_analysis_task_frac = MagmaGeneSetAnalysisTask.create(
        asset_id=base_name + "_drg_cell_analysis_frac",
        magma_gene_analysis_task=gene_analysis_task,
        magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
        gene_set_task=YU_DRG_FRAC_SPECIFICITY_MATRIX,
        set_or_covar="covar",
        model_params=ModelParams(direction_covar="greater", condition_hide=[]),
    )

    cepo_bar_plot = MAGMAPlotGeneSetResult.create(
        gene_set_analysis_task=cell_analysis_task_cepo,
        asset_id=base_name + "_drg_cepo_bar_plot",
        number_of_bars=15,
        label_col="VARIABLE",
    )

    frac_bar_plot = MAGMAPlotGeneSetResult.create(
        gene_set_analysis_task=cell_analysis_task_frac,
        asset_id=base_name + "_drg_frac_bar_plot",
        number_of_bars=15,
        label_col="VARIABLE",
    )
    return YuDRGMAGMA(
        cepo_plot=cepo_bar_plot,
        frac_plot=frac_bar_plot,
    )
