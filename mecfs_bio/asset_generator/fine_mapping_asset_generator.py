"""
Asset generator for fine mapping with SUSIE
"""

import structlog
from attrs import frozen

from mecfs_bio.asset_generator.ukbb_broad_ld_matrix_generator import (
    get_genomic_interval_stem_name,
    get_ld_labels_and_matrix_task_for_genomic_interval_build_37,
    get_optimal_ukbb_ld_interval,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.raw.magma_ensembl_gene_location_reference_data_build_37 import (
    MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_chrom_pos_alleles import (
    HarmonizeGWASWithReferenceViaAlleles,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    ParquetOutFormat,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.build_system.task.pipes.uniquepipe import UniquePipe
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    BroadInstituteFormatLDMatrix,
    SusieRFinemapTask,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    HeatmapOptions,
    RegionSelectDefault,
    SusieStackPlotTask,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)

logger = structlog.get_logger()


@frozen
class BroadFineMapTaskGroup:
    ld_labels_task: Task
    ld_matrix_task: Task
    renamed_ld_labels_task: Task
    harmonized_sumstats_task: Task
    susie_finemap_task: Task
    susie_stackplot_task: Task
    susie_finemap_strict_task: Task
    susie_finemap_strict_plot: Task
    susie_finemap_1_credible_set_task: Task
    susie_finemap_1_credible_set_plot: Task

    def terminal_tasks(self) -> list[Task]:
        return [
            self.susie_finemap_task,
            self.susie_stackplot_task,
            self.susie_finemap_strict_task,
            self.susie_finemap_strict_plot,
            self.susie_finemap_1_credible_set_task,
            self.susie_finemap_1_credible_set_plot,
        ]


def generate_assets_broad_ukbb_fine_map(
    chrom: int,
    pos: int,
    build_37_sumstats_task: Task,
    base_name: str,
    sumstats_pipe: DataProcessingPipe,
    sample_size_or_effect_sample_size: int,
) -> BroadFineMapTaskGroup:
    """
    Asset generator for fine mapping using SUSIE.
    by default, calls SUSIE with a number of different parameter settings to check how the results are affected
    """
    interval = get_optimal_ukbb_ld_interval(
        chrom=chrom,
        pos=pos,
    )
    logger.debug(
        f"To finemap position {pos} on chromosome {chrom}, interval {interval} was selected."
    )
    ld_labels_task, ld_matrix_task = (
        get_ld_labels_and_matrix_task_for_genomic_interval_build_37(
            interval=interval,
        )
    )

    stem = get_genomic_interval_stem_name(interval)
    base_name = base_name + "_" + stem
    ld_labels_task_renamed = PipeDataFrameTask.create(
        source_task=ld_labels_task,
        asset_id=ld_labels_task.asset_id + "_renamed",
        out_format=ParquetOutFormat(),
        pipes=[
            RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
            RenameColPipe(old_name="chromosome", new_name=GWASLAB_CHROM_COL),
            RenameColPipe(old_name="position", new_name=GWASLAB_POS_COL),
            RenameColPipe(
                old_name="allele1",
                new_name=GWASLAB_NON_EFFECT_ALLELE_COL,
                # See: https://github.com/omerwe/polyfun/issues/208#issuecomment-2563832487
            ),
            RenameColPipe(old_name="allele2", new_name=GWASLAB_EFFECT_ALLELE_COL),
        ],
        backend="polars",
    )

    harmonized_sumstats_task = HarmonizeGWASWithReferenceViaAlleles.create(
        asset_id=base_name + "_gwas_harmonized_with_ref",
        gwas_data_task=build_37_sumstats_task,
        reference_task=ld_labels_task_renamed,
        palindrome_strategy="drop",
        gwas_pipe=CompositePipe(
            [
                sumstats_pipe,
                UniquePipe(
                    by=[
                        GWASLAB_CHROM_COL,
                        GWASLAB_POS_COL,
                        GWASLAB_EFFECT_ALLELE_COL,
                        GWASLAB_NON_EFFECT_ALLELE_COL,
                    ],
                    keep="none",
                    order_by=[
                        GWASLAB_CHROM_COL,
                        GWASLAB_POS_COL,
                        GWASLAB_EFFECT_ALLELE_COL,
                        GWASLAB_NON_EFFECT_ALLELE_COL,
                    ],
                ),
            ]
        ),
    )
    susie_finemap_task = SusieRFinemapTask.create(
        asset_id=base_name + "_susie_finemap",
        gwas_data_task=harmonized_sumstats_task,
        ld_labels_task=ld_labels_task_renamed,
        ld_matrix_source=BroadInstituteFormatLDMatrix(ld_matrix_task),
        effective_sample_size=sample_size_or_effect_sample_size,
    )
    # copy_cs=CopyFileFromDirectoryTask.create_result_table(
    #     asset_id=base_name + "_copy_cs_from_directory",
    #     source_directory_task=susie_finemap_task,
    #     path_inside_directory=Path(COMBINED_CS_FILENAME),
    #     extension=".parquet",
    #     read_spec=DataFrameReadSpec(
    #         DataFrameParquetFormat()
    #     )
    #
    # )
    # cs_markdown=ConvertDataFrameToMarkdownTask.create_from_result_table_task(
    #     asset_id=base_name + "_convert_cs_to_markdown",
    #     source_task=copy_cs,
    # )
    susie_stack_plot_task = SusieStackPlotTask.create(
        asset_id=base_name + "_susie_stackplot",
        susie_task=susie_finemap_task,
        gene_info_task=MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
        gene_info_pipe=IdentityPipe(),
        region_mode=RegionSelectDefault(),
        heatmap_options=HeatmapOptions(
            heatmap_bin_options=None, mode="ld2", cmap="plasma"
        ),
    )

    susie_finemap_task_strict = SusieRFinemapTask.create(
        asset_id=base_name + "_susie_finemap_strict_threshold",
        gwas_data_task=harmonized_sumstats_task,
        ld_labels_task=ld_labels_task_renamed,
        ld_matrix_source=BroadInstituteFormatLDMatrix(ld_matrix_task),
        effective_sample_size=sample_size_or_effect_sample_size,
        z_score_filtering_threshold=1.0,
    )
    strict_plot = SusieStackPlotTask.create(
        asset_id=base_name + "_susie_stackplot_strict",
        susie_task=susie_finemap_task_strict,
        gene_info_task=MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
        gene_info_pipe=IdentityPipe(),
        region_mode=RegionSelectDefault(),
        heatmap_options=HeatmapOptions(
            heatmap_bin_options=None, mode="ld2", cmap="plasma"
        ),
    )

    susie_finemap_task_1_credible_set = SusieRFinemapTask.create(
        asset_id=base_name + "_susie_finemap_1_credible_set",
        gwas_data_task=harmonized_sumstats_task,
        ld_labels_task=ld_labels_task_renamed,
        ld_matrix_source=BroadInstituteFormatLDMatrix(ld_matrix_task),
        effective_sample_size=sample_size_or_effect_sample_size,
        max_credible_sets=1,
        # z_score_filtering_threshold=1.0
    )
    susie_plot_1_credible_set = SusieStackPlotTask.create(
        asset_id=base_name + "_susie_stackplot_1_credible_set",
        susie_task=susie_finemap_task_1_credible_set,
        gene_info_task=MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
        gene_info_pipe=IdentityPipe(),
        region_mode=RegionSelectDefault(),
        heatmap_options=HeatmapOptions(
            heatmap_bin_options=None, mode="ld2", cmap="plasma"
        ),
    )

    return BroadFineMapTaskGroup(
        ld_labels_task=ld_labels_task,
        ld_matrix_task=ld_matrix_task,
        renamed_ld_labels_task=ld_labels_task_renamed,
        harmonized_sumstats_task=harmonized_sumstats_task,
        susie_finemap_task=susie_finemap_task,
        susie_stackplot_task=susie_stack_plot_task,
        susie_finemap_strict_task=susie_finemap_task_strict,
        susie_finemap_strict_plot=strict_plot,
        susie_finemap_1_credible_set_task=susie_finemap_task_1_credible_set,
        susie_finemap_1_credible_set_plot=susie_plot_1_credible_set,
    )
