from typing import Mapping, Sequence

from attrs import frozen

from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.fdr_multiple_testing_table_task import (
    Method,
    MultipleTestingTableTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_cell_analysis_by_ldsc import (
    CellAnalysisByLDSCTask,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.sldsc_scatter_plot_task import SLDSCScatterPlotTask


@frozen
class PartitionedLDScoresRecord:
    entry_name: str
    ref_ld_chr_cts_task: Task
    ref_ld_chr_cts_filename: str
    cell_or_tissue_labels_task: Task | None


@frozen
class CellAnalysisTaskGroup:
    cell_analysis_task: CellAnalysisByLDSCTask
    multiple_testing_task: Task
    add_categories_task: Task | None
    plot_task: Task | None

    def terminal_task(self) -> Task:
        if self.plot_task is not None:
            return self.plot_task
        if self.add_categories_task is not None:
            return self.add_categories_task
        return self.multiple_testing_task


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
        return [item.terminal_task() for item in self.partitioned_tasks.values()]

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
            )
            multiple_testing_task = (
                MultipleTestingTableTask.create_from_result_table_task(
                    asset_id=base_name
                    + "_"
                    + entry.entry_name
                    + "_multiple_testing_correction",
                    alpha=multiple_testing_alpha,
                    p_value_column="Coefficient_P_value",
                    source_task=cell_task,
                    apply_filter=False,
                    method=multiple_testing_method,
                )
            )
            if entry.cell_or_tissue_labels_task is not None:
                add_labels_task = JoinDataFramesTask.create_from_result_df(
                    asset_id=base_name + "_" + entry.entry_name + "_add_labels",
                    result_df_task=multiple_testing_task,
                    reference_df_task=entry.cell_or_tissue_labels_task,
                    how="left",
                    left_on=["Name"],
                    right_on=["Tissue_Or_Cell"],
                )
                plot_task = SLDSCScatterPlotTask.create(
                    asset_id=base_name + "_" + entry.entry_name + "_sldsc_plot",
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
            )
        return cls(cell_analysis_task_groups)
