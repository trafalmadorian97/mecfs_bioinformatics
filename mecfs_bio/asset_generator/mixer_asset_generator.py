"""
Asset generator to use the MiXeR Gaussian mixture to model the genetic architecture of a trait.
"""

from pathlib import PurePath
from typing import Mapping, Sequence

from attrs import frozen

from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.convert_dataframe_to_markdown_task import (
    ConvertDataFrameToMarkdownTask,
)
from mecfs_bio.build_system.task.copy_file_from_directory_task import (
    CopyFileFromDirectoryTask,
)
from mecfs_bio.build_system.task.mixer.mixer_task import (
    MixerDataSource,
    MixerTask,
    UnivariateMode,
)
from mecfs_bio.build_system.task.mixer.mixer_univariate_combine import (
    MixerRunSource,
    MixerUnivariateCombine,
)
from mecfs_bio.build_system.task.mixer.mixer_univariate_results import (
    TEST_OUTPUT_PREFIX,
    MixerUnivariateSummarizeResultsTask,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.drop_col_pipe import DropColPipe
from mecfs_bio.build_system.task.pipes.heritability_conversion_pipe import (
    HeritabilityConversionPipe,
)
from mecfs_bio.build_system.task.pipes.rename_col_by_position_pipe import (
    RenameColByPositionPipe,
)
from mecfs_bio.build_system.task.pipes.transpose_pipe import TransposePipe


@frozen
class UnivariateMixerTasks:
    run_tasks: Mapping[int, Task]
    combine_task: Task
    results_task: Task
    power_plot_task: Task
    qq_plot_task: Task
    qq_bin_plot_task: Task
    result_markdown_table_task: Task

    def terminal_tasks(self) -> list[Task]:
        return [
            self.results_task,
            self.power_plot_task,
            self.qq_plot_task,
            self.qq_bin_plot_task,
            self.result_markdown_table_task,
        ]


def univariate_mixer_asset_generator(
    base_name: str,
    name_in_plot: str,
    trait_1_source: MixerDataSource,
    reference_data_directory_task: Task,
    reps: Sequence[int] = tuple(range(1, 21)),
    threads: int = 4,
):
    """
    Asset generator to apply univariate MiXeR to GWAS summary statistics.

    See:
    Holland, Dominic, et al. "Beyond SNP heritability: Polygenicity and discoverability of phenotypes
    estimated with a univariate Gaussian mixture model." PLoS Genetics 16.5 (2020): e1008612.
    """
    tasks = {}
    for rep in reps:
        tasks[rep] = MixerTask.create(
            asset_id=base_name + f"_standard_mixer_rep_{rep}",
            trait_1_source=trait_1_source,
            mixer_mode=UnivariateMode(),
            ref_data_directory_task=reference_data_directory_task,
            reps_to_perform=[rep],
            threads=threads,
        )
    combine_task = MixerUnivariateCombine.create(
        asset_id=base_name + "_univariate_mixer_combine",
        mixer_source_runs=[
            MixerRunSource(task=task, rep=num) for num, task in tasks.items()
        ],
        trait_name=name_in_plot,
    )
    result_task = MixerUnivariateSummarizeResultsTask.create(
        asset_id=base_name + "_univariate_mixer_results_dir",
        combine_task=combine_task,
        trait_name=name_in_plot,
    )
    power_plot_task = CopyFileFromDirectoryTask.create_from_result_plot(
        asset_id=base_name + "_univariate_mixer_power_plot",
        source_directory_task=result_task,
        path_inside_directory=PurePath(TEST_OUTPUT_PREFIX + ".power.png"),
        extension=".png",
        subdir=PurePath("analysis/mixer_plots"),
    )
    qq_plot_task = CopyFileFromDirectoryTask.create_from_result_plot(
        asset_id=base_name + "_univariate_mixer_qq_plot",
        source_directory_task=result_task,
        path_inside_directory=PurePath(TEST_OUTPUT_PREFIX + ".qq.png"),
        extension=".png",
        subdir=PurePath("analysis/mixer_plots"),
    )
    qq_bin_plot_task = CopyFileFromDirectoryTask.create_from_result_plot(
        asset_id=base_name + "_univariate_mixer_qqbin_plot",
        source_directory_task=result_task,
        path_inside_directory=PurePath(TEST_OUTPUT_PREFIX + ".qqbin.png"),
        extension=".png",
        subdir=PurePath("analysis/mixer_plots"),
    )
    result_table_task = CopyFileFromDirectoryTask.create_result_table(
        asset_id=base_name + "_univariate_mixer_results_table",
        source_directory_task=result_task,
        path_inside_directory=PurePath(TEST_OUTPUT_PREFIX + ".csv"),
        extension=".csv",
        read_spec=DataFrameReadSpec(DataFrameTextFormat(separator="\t")),
    )

    result_table_as_markdown_task = (
        ConvertDataFrameToMarkdownTask.create_from_result_table_task(
            source_task=result_table_task,
            asset_id=base_name + "_univariate_mixer_results_table_as_markdown",
            pipe=CompositePipe(
                [
                    DropColPipe(["fname"]),
                    HeritabilityConversionPipe(
                        observed_heritability_column="h2 (mean)",
                        liability_heritability_column="h2 (liability, mean)",
                        sample_info=trait_1_source.sample_info,
                        assume_even_sample=True,
                    ),
                    TransposePipe(),
                    RenameColByPositionPipe(col_position=0, col_new_name="Parameter"),
                    RenameColByPositionPipe(col_position=1, col_new_name="Value"),
                ]
            ),
        )
    )
    return UnivariateMixerTasks(
        tasks,
        combine_task=combine_task,
        results_task=result_task,
        power_plot_task=power_plot_task,
        qq_plot_task=qq_plot_task,
        qq_bin_plot_task=qq_bin_plot_task,
        result_markdown_table_task=result_table_as_markdown_task,
    )
