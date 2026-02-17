"""
Asset generator for fine mapping with SUSIE
"""

from pathlib import Path, PurePath

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
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.convert_dataframe_to_markdown_task import (
    ConvertDataFrameToMarkdownTask,
)
from mecfs_bio.build_system.task.copy_file_from_directory_task import (
    CopyFileFromDirectoryTask,
)
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_chrom_pos_alleles import (
    ChromRange,
    HarmonizeGWASWithReferenceViaAlleles,
)
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_rsid import (
    PalindromeStrategy,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    ParquetOutFormat,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.concat_str_pipe import ConcatStrPipe
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.filter_rows_by_min_in_col import (
    FilterRowsByMinInCol,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.build_system.task.pipes.uniquepipe import UniquePipe
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    COMBINED_CS_FILENAME,
    PIP_COLUMN,
    BroadInstituteFormatLDMatrix,
    SusieRFinemapTask,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    HeatmapOptions,
    RegionSelectDefault,
    SusieStackPlotTask,
)
from mecfs_bio.build_system.task.upset_plot_task import DirSetSource, UpSetPlotTask
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
    susie_finemap_2_credible_set_task: Task
    susie_finemap_2_credible_set_plot: Task
    markdown_table_tasks: list[Task]
    upset_plot_task: Task
    upset_plot_task_pip001: Task

    def terminal_tasks(self) -> list[Task]:
        return [
            self.susie_finemap_task,
            self.susie_stackplot_task,
            self.susie_finemap_strict_task,
            self.susie_finemap_strict_plot,
            self.susie_finemap_1_credible_set_task,
            self.susie_finemap_1_credible_set_plot,
            self.susie_finemap_2_credible_set_task,
            self.susie_finemap_2_credible_set_plot,
            self.upset_plot_task,
            # self.upset_plot_task_pip005,
            # self.upset_plot_task_pip0025,
            self.upset_plot_task_pip001,
        ] + self.markdown_table_tasks


def generate_assets_broad_ukbb_fine_map(
    chrom: int,
    pos: int,
    build_37_sumstats_task: Task,
    base_name: str,
    sumstats_pipe: DataProcessingPipe,
    sample_size_or_effect_sample_size: int,
    chrom_range: ChromRange | None = None,
    palindrome_strategy: PalindromeStrategy = "drop",
) -> BroadFineMapTaskGroup:
    """
    Asset generator for fine mapping using SUSIE.
    by default, calls SUSIE with a number of different parameter settings to check how the results are affected
    """

    interval = get_optimal_ukbb_ld_interval(
        chrom=chrom,
        pos=pos,
    )
    if chrom_range is not None:
        assert chrom == chrom_range.chrom
        assert pos >= chrom_range.start
        assert pos <= chrom_range.end
        assert interval.start <= chrom_range.start
        assert chrom_range.end <= interval.end
        base_name = (
            base_name + f"chr{chrom_range.chrom}_{chrom_range.start}_{chrom_range.end}"
        )
    else:
        stem = get_genomic_interval_stem_name(interval)
        base_name = base_name + "_" + stem
    if palindrome_strategy != "drop":
        base_name = base_name + "_palindromes_" + palindrome_strategy

    logger.debug(
        f"To finemap position {pos} on chromosome {chrom}, interval {interval} was selected."
    )
    ld_labels_task, ld_matrix_task = (
        get_ld_labels_and_matrix_task_for_genomic_interval_build_37(
            interval=interval,
        )
    )

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
        palindrome_strategy=palindrome_strategy,
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
        chrom_range_filter=chrom_range,
    )
    markdown_table_tasks = []
    susie_finemap_task = SusieRFinemapTask.create(
        asset_id=base_name + "_susie_finemap",
        gwas_data_task=harmonized_sumstats_task,
        ld_labels_task=ld_labels_task_renamed,
        ld_matrix_source=BroadInstituteFormatLDMatrix(ld_matrix_task),
        effective_sample_size=sample_size_or_effect_sample_size,
    )
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
    markdown_table_tasks.append(
        markdown_cs_table_task(
            susie_finemap_task=susie_finemap_task, base_name=base_name + "_susie_base"
        )
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

    markdown_table_tasks.append(
        markdown_cs_table_task(
            susie_finemap_task=susie_finemap_task_strict,
            base_name=base_name + "_susie_strict",
        )
    )

    susie_finemap_task_1_credible_set = SusieRFinemapTask.create(
        asset_id=base_name + "_susie_finemap_1_credible_set",
        gwas_data_task=harmonized_sumstats_task,
        ld_labels_task=ld_labels_task_renamed,
        ld_matrix_source=BroadInstituteFormatLDMatrix(ld_matrix_task),
        effective_sample_size=sample_size_or_effect_sample_size,
        max_credible_sets=1,
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

    markdown_table_tasks.append(
        markdown_cs_table_task(
            susie_finemap_task=susie_finemap_task_1_credible_set,
            base_name=base_name + "_susie_1",
        )
    )

    susie_finemap_task_2_credible_set = SusieRFinemapTask.create(
        asset_id=base_name + "_susie_finemap_2_credible_set",
        gwas_data_task=harmonized_sumstats_task,
        ld_labels_task=ld_labels_task_renamed,
        ld_matrix_source=BroadInstituteFormatLDMatrix(ld_matrix_task),
        effective_sample_size=sample_size_or_effect_sample_size,
        max_credible_sets=1,
    )
    susie_plot_2_credible_set = SusieStackPlotTask.create(
        asset_id=base_name + "_susie_stackplot_2_credible_set",
        susie_task=susie_finemap_task_2_credible_set,
        gene_info_task=MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
        gene_info_pipe=IdentityPipe(),
        region_mode=RegionSelectDefault(),
        heatmap_options=HeatmapOptions(
            heatmap_bin_options=None, mode="ld2", cmap="plasma"
        ),
    )

    markdown_table_tasks.append(
        markdown_cs_table_task(
            susie_finemap_task=susie_finemap_task_2_credible_set,
            base_name=base_name + "_susie_2",
        )
    )
    variant_id = "__variant_id"
    id_variant_pipe = ConcatStrPipe(
        target_cols=[
            GWASLAB_CHROM_COL,
            GWASLAB_POS_COL,
            GWASLAB_EFFECT_ALLELE_COL,
            GWASLAB_NON_EFFECT_ALLELE_COL,
        ],
        sep="__",
        new_col_name=variant_id,
    )

    filtered_id_variant_pipe_005 = CompositePipe(
        [id_variant_pipe, FilterRowsByMinInCol(min_value=0.05, col=PIP_COLUMN)]
    )

    filtered_id_variant_pipe_0025 = CompositePipe(
        [id_variant_pipe, FilterRowsByMinInCol(min_value=0.025, col=PIP_COLUMN)]
    )

    filtered_id_variant_pipe_001 = CompositePipe(
        [id_variant_pipe, FilterRowsByMinInCol(min_value=0.01, col=PIP_COLUMN)]
    )

    upset_plot = UpSetPlotTask.create(
        asset_id=base_name + "_upset_plot",
        set_sources=[
            DirSetSource(
                name="L=10",
                task=susie_finemap_task,
                file_in_dir=PurePath(COMBINED_CS_FILENAME),
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
                pipe=id_variant_pipe,
                col_name=variant_id,
            ),
            DirSetSource(
                name="L=10, Strict",
                task=susie_finemap_task_strict,
                file_in_dir=PurePath(COMBINED_CS_FILENAME),
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
                pipe=id_variant_pipe,
                col_name=variant_id,
            ),
            DirSetSource(
                name="L=2",
                task=susie_finemap_task_2_credible_set,
                file_in_dir=PurePath(COMBINED_CS_FILENAME),
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
                pipe=id_variant_pipe,
                col_name=variant_id,
            ),
            DirSetSource(
                name="L=1",
                task=susie_finemap_task_1_credible_set,
                file_in_dir=PurePath(COMBINED_CS_FILENAME),
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
                pipe=id_variant_pipe,
                col_name=variant_id,
            ),
        ],
    )

    upset_plot_pip001 = UpSetPlotTask.create(
        asset_id=base_name + "_upset_plot_pip001",
        set_sources=[
            DirSetSource(
                name="L=10",
                task=susie_finemap_task,
                file_in_dir=PurePath(COMBINED_CS_FILENAME),
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
                pipe=filtered_id_variant_pipe_001,
                col_name=variant_id,
            ),
            DirSetSource(
                name="L=10, Strict",
                task=susie_finemap_task_strict,
                file_in_dir=PurePath(COMBINED_CS_FILENAME),
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
                pipe=filtered_id_variant_pipe_001,
                col_name=variant_id,
            ),
            DirSetSource(
                name="L=2",
                task=susie_finemap_task_2_credible_set,
                file_in_dir=PurePath(COMBINED_CS_FILENAME),
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
                pipe=filtered_id_variant_pipe_001,
                col_name=variant_id,
            ),
            DirSetSource(
                name="L=1",
                task=susie_finemap_task_1_credible_set,
                file_in_dir=PurePath(COMBINED_CS_FILENAME),
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
                pipe=filtered_id_variant_pipe_001,
                col_name=variant_id,
            ),
        ],
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
        markdown_table_tasks=markdown_table_tasks,
        susie_finemap_2_credible_set_task=susie_finemap_task_2_credible_set,
        susie_finemap_2_credible_set_plot=susie_plot_2_credible_set,
        upset_plot_task=upset_plot,
        upset_plot_task_pip001=upset_plot_pip001,
    )


def markdown_cs_table_task(
    susie_finemap_task: SusieRFinemapTask, base_name: str
) -> Task:
    copy_cs = CopyFileFromDirectoryTask.create_result_table(
        asset_id=base_name + "_copy_cs_from_directory",
        source_directory_task=susie_finemap_task,
        path_inside_directory=Path(COMBINED_CS_FILENAME),
        extension=".parquet",
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    )

    cs_markdown = ConvertDataFrameToMarkdownTask.create_from_result_table_task(
        asset_id=base_name + "_convert_cs_to_markdown",
        source_task=copy_cs,
    )
    return cs_markdown
