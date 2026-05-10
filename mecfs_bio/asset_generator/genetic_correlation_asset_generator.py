"""
Task generator that, given a list of tasks that generate gwaslab sumstats,
creates one task for each pair of sumstats to compute their genetic correlation.
The advantage of this approach is that each correlation is cached separately.
"""

import itertools
from typing import Mapping, Sequence

import narwhals.dtypes
from attrs import frozen
from frozendict import frozendict

from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genome_phase_3_v1_extracted import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.concat_frames_task import ConcatFramesTask
from mecfs_bio.build_system.task.convert_dataframe_to_markdown_task import (
    ConvertDataFrameToMarkdownTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    GeneticCorrelationByCTLDSCTask,
    SumstatsSource,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_snp_heritability_by_ldsc_task import (
    SNPHeritabilityByLDSCTask,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import CSVOutFormat
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.drop_col_pipe import DropColPipe
from mecfs_bio.build_system.task.pipes.format_numbers_pipe import FormatFloatNumbersPipe
from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild


@frozen
class GeneticCorrTasks:
    corr_tasks: Mapping[str, Task]
    aggregation_task: Task
    heritability_tasks: Sequence[Task]
    heritability_aggregation_task: Task
    aggregation_markdown_task: Task

    def terminal_tasks(self) -> list[Task]:
        return [
            self.aggregation_task,
            self.heritability_aggregation_task,
            self.aggregation_markdown_task,
        ]


def genetic_corr_by_ct_ldsc_asset_generator(
    base_name: str,
    sources: Sequence[SumstatsSource],
    ld_ref_task: Task = THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
    ld_file_filename_pattern: str = "/LDscore.@",
    build: GenomeBuild = "19",
) -> GeneticCorrTasks:
    """
    Task generator that, given a list of tasks that generate gwaslab sumstats, creates one task for each pair of sumstats to compute their genetic correlation via CT-LDSC.
    The advantage of this approach is that each correlation is cached separately.
    """
    tasks = {}
    for source_1, source_2 in itertools.combinations(sources, 2):
        sorted_sources = sorted([source_1, source_2], key=lambda x: x.alias.lower())
        task, task_name = _create_task_for_source_pair(
            sorted_sources=sorted_sources,
            base_name=base_name,
            ld_ref_task=ld_ref_task,
            ld_file_filename_pattern=ld_file_filename_pattern,
            build=build,
        )
        tasks[task_name] = task

    # Create the aggregation task to combine all the covariances
    aggregation = ConcatFramesTask.create(
        asset_id=base_name + "_ct_ldsc_corr_aggregation",
        frames_tasks=list(tasks.values()),
        out_format=CSVOutFormat(","),
    )

    aggregation_md = ConvertDataFrameToMarkdownTask.create_from_result_table_task(
        asset_id=base_name + "_ct_ldsc_corr_agg_markdown",
        source_task=aggregation,
        pipe=CompositePipe(
            [DropColPipe(["h2_liab", "h2_liab_se", "h2_int", "h2_int_se"])]
            + [
                # TransposePipe(),
                # RenameColByPositionPipe(
                #     col_position=0, col_new_name="Parameter"
                # ),
                # RenameColByPositionPipe(col_position=1, col_new_name="Value"),
                # SelectColPipe(["Parameter", "Value"]),
                FormatFloatNumbersPipe(col=item, format_str=".4g")
                for item in ["rg", "se", "z", "p", "gcov_int", "gcov_int_se"]
            ]
        ),
    )

    heritability_tasks = [
        SNPHeritabilityByLDSCTask.create(
            asset_id=base_name + "_heritability_" + source.alias,
            source_sumstats_task=source.task,
            phenotype_info=source.sample_info,
            pipe=source.pipe,
            ld_ref_task=ld_ref_task,
            ld_file_filename_pattern=ld_file_filename_pattern,
            build=build,
        )
        for source in sources
    ]
    heritability_agg = ConcatFramesTask.create(
        asset_id=base_name + "_heritability_aggregation",
        frames_tasks=heritability_tasks,
        out_format=CSVOutFormat(","),
        override_trait="multi_trait",
        override_project="heritability",
        column_type_override={
            "Ratio": narwhals.dtypes.String(),
            "Ratio_se": narwhals.dtypes.String(),
        },
        frames_pipes=[
            SetColToConstantPipe(col_name="p", constant=source.alias)
            for source in sources
        ],
    )

    return GeneticCorrTasks(
        corr_tasks=frozendict(tasks),
        aggregation_task=aggregation,
        heritability_tasks=heritability_tasks,
        heritability_aggregation_task=heritability_agg,
        aggregation_markdown_task=aggregation_md,
    )


def _create_task_for_source_pair(
    sorted_sources: Sequence[SumstatsSource],
    base_name: str,
    ld_ref_task: Task,
    build: GenomeBuild,
    ld_file_filename_pattern: str,
) -> tuple[Task, str]:
    assert len(sorted_sources) == 2
    task_name = (
        base_name
        + "_ct_ldsc_corr_"
        + sorted_sources[0].alias
        + "__"
        + sorted_sources[1].alias
    ).lower()
    task = GeneticCorrelationByCTLDSCTask.create(
        asset_id=task_name,
        sources=sorted_sources,
        ld_ref_task=ld_ref_task,
        build=build,
        ld_file_filename_pattern=ld_file_filename_pattern,
    )
    return task, task_name
