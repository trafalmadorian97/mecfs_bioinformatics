from typing import Mapping, Self, Sequence

from attrs import frozen

from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.convert_dataframe_to_markdown_task import (
    ConvertDataFrameToMarkdownTask,
)
from mecfs_bio.build_system.task.dataframe_output import (
    ParquetOutFormat,
    ParquetWriteOptions,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_cell_analysis_by_sldsc import (
    CellAnalysisByLDSCTask,
)
from mecfs_bio.build_system.task.gwaslab.sldsc_scatter_plot_task import (
    SLDSCScatterPlotTask,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.multiple_testing_table_task import (
    Method,
    MultipleTestingTableTask,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import PipeDataFrameTask
from mecfs_bio.build_system.task.pipes.cast_to_int32pipe import CastIntsToInt32Pipe
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.drop_col_pipe import DropColPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.constants.sldsc_constants import (
    SLDSC_CTS_COEF_SE_COL,
    SLDSC_CTS_P_COL,
)
from mecfs_bio.util.type_related.unwrap import unwrap

# Appended to every asset id generated for the vectorized S-LDSC Task, so that its results
# coexist with those of CellAnalysisByLDSCTask.
VEC_ID_SUFFIX = "_vec"

# Encoding for tables shown on the documentation site via the data_table macro.
# Int32 rather than Int64 because parquet INT64 decodes to JavaScript BigInt.
_DOCS_PARQUET_FORMAT = ParquetOutFormat(
    write_options=ParquetWriteOptions(
        compression="zstd",
        compression_level=22,
        byte_stream_split_floats=True,
    )
)


@frozen(slots=True)
class CellOrTissueLabelRecord:
    cell_or_tissue_label_task: Task
    pipe_left: DataProcessingPipe = IdentityPipe()
    pipe_right: DataProcessingPipe = IdentityPipe()
    left_join_on: str = "Name"
    right_join_on: str = "Tissue_Or_Cell"


@frozen(slots=True)
class PartitionedLDScoresRecord:
    entry_name: str
    ref_ld_chr_cts_task: Task
    ref_ld_chr_cts_filename: str
    cell_or_tissue_labels_task: CellOrTissueLabelRecord | None


@frozen(slots=True)
class CellAnalysisTaskGroup:
    cell_analysis_task: Task
    multiple_testing_task: Task
    multiple_testing_task_markdown: Task
    multiple_testing_task_parquet: Task
    add_categories_task: Task | None
    plot_task: Task | None

    @property
    def plot_task_unwrap(self) -> Task:
        return unwrap(self.plot_task)

    @property
    def add_categories_task_unwrap(self) -> Task:
        return unwrap(self.add_categories_task)

    def terminal_tasks(self) -> list[Task]:
        result = [
            self.multiple_testing_task_markdown,
            self.multiple_testing_task_parquet,
        ]
        if self.plot_task is not None:
            result.append(self.plot_task)
        return result


@frozen(slots=True)
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
    ) -> Self:
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
            cell_analysis_task_groups[entry.entry_name] = (
                build_cell_analysis_task_group(
                    base_name=base_name,
                    entry_name=entry.entry_name,
                    cell_task=cell_task,
                    labels=entry.cell_or_tissue_labels_task,
                    multiple_testing_alpha=multiple_testing_alpha,
                    multiple_testing_method=multiple_testing_method,
                    id_suffix="",
                )
            )
        return cls(cell_analysis_task_groups)


def build_cell_analysis_task_group(
    base_name: str,
    entry_name: str,
    cell_task: Task,
    labels: CellOrTissueLabelRecord | None,
    multiple_testing_alpha: float,
    multiple_testing_method: Method,
    id_suffix: str,
) -> CellAnalysisTaskGroup:
    """
    Downstream tasks for one cell-type-specific S-LDSC result table: multiple-testing
    correction, markdown and parquet presentation tables, and (when labels are given) a label
    join and scatter plot. id_suffix is appended to every generated asset id.
    """
    prefix = base_name + "_" + entry_name
    multiple_testing_task = MultipleTestingTableTask.create_from_result_table_task(
        asset_id=prefix + "_s_ldsc_multiple_testing_correction" + id_suffix,
        alpha=multiple_testing_alpha,
        p_value_column=SLDSC_CTS_P_COL,
        source_task=cell_task,
        apply_filter=False,
        method=multiple_testing_method,
    )
    presentation_drop_cols = DropColPipe(["_Corrected P Value_", SLDSC_CTS_COEF_SE_COL])
    multiple_testing_markdown = (
        ConvertDataFrameToMarkdownTask.create_from_result_table_task(
            source_task=multiple_testing_task,
            asset_id=prefix + "_s_ldsc_cell_analysis_md_table" + id_suffix,
            pipe=presentation_drop_cols,
        )
    )
    multiple_testing_parquet = PipeDataFrameTask.create(
        source_task=multiple_testing_task,
        asset_id=prefix + "_s_ldsc_cell_analysis_parquet_table" + id_suffix,
        pipes=[presentation_drop_cols, CastIntsToInt32Pipe()],
        out_format=_DOCS_PARQUET_FORMAT,
        backend="polars",
    )
    if labels is not None:
        add_labels_task = JoinDataFramesTask.create_from_result_df(
            asset_id=prefix + "_add_labels" + id_suffix,
            result_df_task=multiple_testing_task,
            reference_df_task=labels.cell_or_tissue_label_task,
            how="left",
            left_on=[labels.left_join_on],
            right_on=[labels.right_join_on],
            df_1_pipe=labels.pipe_left,
            df_2_pipe=labels.pipe_right,
        )
        plot_task = SLDSCScatterPlotTask.create(
            asset_id=prefix + "_cell_analysis_s_ldsc_plot" + id_suffix,
            source_task=add_labels_task,
        )
    else:
        add_labels_task = None
        plot_task = None
    return CellAnalysisTaskGroup(
        cell_analysis_task=cell_task,
        multiple_testing_task=multiple_testing_task,
        add_categories_task=add_labels_task,
        plot_task=plot_task,
        multiple_testing_task_markdown=multiple_testing_markdown,
        multiple_testing_task_parquet=multiple_testing_parquet,
    )
