from typing import Mapping, Sequence

from attrs import frozen

from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.convert_dataframe_to_markdown_task import (
    ConvertDataFrameToMarkdownTask,
)
from mecfs_bio.build_system.task.fdr_multiple_testing_table_task import (
    Method,
    MultipleTestingTableTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_cell_analysis_by_sldsc import (
    CellAnalysisByLDSCTask,
)
from mecfs_bio.build_system.task.gwaslab.sldsc_scatter_plot_task import (
    SLDSCScatterPlotTask,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.drop_col_pipe import DropColPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe


@frozen
class CellOrTissueLabelRecord:
    cell_or_tissue_label_task: Task
    pipe_left: DataProcessingPipe = IdentityPipe()
    pipe_right: DataProcessingPipe = IdentityPipe()
    left_join_on: str = "Name"
    right_join_on: str = "Tissue_Or_Cell"


@frozen
class PartitionedLDScoresRecord:
    entry_name: str
    ref_ld_chr_cts_task: Task
    ref_ld_chr_cts_filename: str
    cell_or_tissue_labels_task: CellOrTissueLabelRecord | None


@frozen
class CellAnalysisTaskGroup:
    cell_analysis_task: CellAnalysisByLDSCTask
    multiple_testing_task: Task
    multiple_testing_task_markdown: Task
    add_categories_task: Task | None
    plot_task: Task | None

    def terminal_tasks(self) -> list[Task]:
        result = [self.multiple_testing_task_markdown]
        if self.plot_task is not None:
            result.append(self.plot_task)
        return result


@frozen
class SLDSCTaskGenerator:
    """
       Task generator to apply cell- and tissue-type specific Stratified Linkage Disequilibrium Score Regression
       to dataset as described in

    "Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types." Nature genetics 50.4 (2018): 621-629.

       Assumption is that a variety of different partitioning datasets will be applied to a single set of GWAS data.

    """

    partitioned_tasks: Mapping[str, CellAnalysisTaskGroup]

    def get_terminal_tasks(self) -> Sequence[Task]:
        result = []
        for item in self.partitioned_tasks.values():
            result.extend(item.terminal_tasks())
        return result

    @classmethod
    def create(
        cls,
        base_name: str,
        source_sumstats_task: Task,
        ref_ld_chr_task: Task,
        ref_ld_chr_inner_dirname: str,
        w_ld_chr_task: Task,
        w_ld_chr_inner_dirname: str,
        partitioned_entries: list[PartitionedLDScoresRecord],
        multiple_testing_alpha: float,
        multiple_testing_method: Method,
        pre_pipe: DataProcessingPipe = IdentityPipe(),
    ):
        cell_analysis_task_groups = {}
        for entry in partitioned_entries:
            cell_task = CellAnalysisByLDSCTask.create(
                asset_id=base_name
                + "_"
                + entry.entry_name
                + "_cell_analysis_by_s_ldsc",
                source_sumstats_task=source_sumstats_task,
                ref_ld_chr_cts_task=entry.ref_ld_chr_cts_task,
                ref_ld_chr_cts_filename=entry.ref_ld_chr_cts_filename,
                ref_ld_chr_task=ref_ld_chr_task,
                ref_ld_chr_inner_dirname=ref_ld_chr_inner_dirname,
                w_ld_chr_task=w_ld_chr_task,
                w_ld_chr_inner_dirname=w_ld_chr_inner_dirname,
                pre_pipe=pre_pipe,
            )
            multiple_testing_task = (
                MultipleTestingTableTask.create_from_result_table_task(
                    asset_id=base_name
                    + "_"
                    + entry.entry_name
                    + "_s_ldsc_multiple_testing_correction",
                    alpha=multiple_testing_alpha,
                    p_value_column="Coefficient_P_value",
                    source_task=cell_task,
                    apply_filter=False,
                    method=multiple_testing_method,
                )
            )
            multiple_testing_markdown = (
                ConvertDataFrameToMarkdownTask.create_from_result_table_task(
                    source_task=multiple_testing_task,
                    asset_id=base_name
                    + "_"
                    + entry.entry_name
                    + "_s_ldsc_cell_analysis_md_table",
                    pipe=DropColPipe(["_Corrected P Value_", "Coefficient_std_error"]),
                )
            )
            if entry.cell_or_tissue_labels_task is not None:
                add_labels_task = JoinDataFramesTask.create_from_result_df(
                    asset_id=base_name + "_" + entry.entry_name + "_add_labels",
                    result_df_task=multiple_testing_task,
                    reference_df_task=entry.cell_or_tissue_labels_task.cell_or_tissue_label_task,
                    how="left",
                    left_on=[entry.cell_or_tissue_labels_task.left_join_on],
                    right_on=[entry.cell_or_tissue_labels_task.right_join_on],
                    df_1_pipe=entry.cell_or_tissue_labels_task.pipe_left,
                    df_2_pipe=entry.cell_or_tissue_labels_task.pipe_right,
                )
                plot_task = SLDSCScatterPlotTask.create(
                    asset_id=base_name
                    + "_"
                    + entry.entry_name
                    + "_cell_analysis_s_ldsc_plot",
                    source_task=add_labels_task,
                )
            else:
                add_labels_task = None
                plot_task = None
            cell_analysis_task_groups[entry.entry_name] = CellAnalysisTaskGroup(
                cell_analysis_task=cell_task,
                multiple_testing_task=multiple_testing_task,
                add_categories_task=add_labels_task,
                plot_task=plot_task,
                multiple_testing_task_markdown=multiple_testing_markdown,
            )
        return cls(cell_analysis_task_groups)
