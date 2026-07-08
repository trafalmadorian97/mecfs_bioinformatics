"""
Task generator wiring the POPs pipeline: munge gene features, run POPs against a
MAGMA gene-analysis result, and extract the per-gene predictions into a standalone
result table.
"""

from pathlib import PurePath

from attrs import frozen

from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameWhiteSpaceSepTextFormat,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.copy_file_from_directory_task import (
    CopyFileFromDirectoryTask,
)
from mecfs_bio.build_system.task.pops.pops_munge_task import (
    PopsMungeFeatureDirectoryTask,
)
from mecfs_bio.build_system.task.pops.pops_run_task import PopsRunTask
from mecfs_bio.build_system.task.pops.pops_utils import (
    POPS_OUTPUT_STEM_NAME,
    POPS_PREDS_SUFFIX,
)


@frozen
class PopsTaskGenerator:
    """The Tasks required to run POPs on top of a MAGMA gene analysis."""

    munge_task: Task
    pops_run_task: Task
    preds_result_table_task: Task

    def terminal_tasks(self) -> list[Task]:
        return [self.preds_result_table_task]

    @classmethod
    def create(
        cls,
        base_name: str,
        pops_source_task: Task,
        raw_features_task: Task,
        magma_gene_analysis_task: Task,
        max_cols: int = 500,
        pops_extra_args: tuple[str, ...] = (),
    ) -> "PopsTaskGenerator":
        munge_task = PopsMungeFeatureDirectoryTask.create(
            asset_id=base_name + "_pops_munge_features",
            pops_source_task=pops_source_task,
            raw_features_task=raw_features_task,
            max_cols=max_cols,
        )
        pops_run_task = PopsRunTask.create(
            asset_id=base_name + "_pops_run",
            pops_source_task=pops_source_task,
            munged_features_task=munge_task,
            magma_gene_analysis_task=magma_gene_analysis_task,
            extra_args=pops_extra_args,
        )
        preds_result_table_task = CopyFileFromDirectoryTask.create_result_table(
            asset_id=base_name + "_pops_preds",
            source_directory_task=pops_run_task,
            path_inside_directory=PurePath(POPS_OUTPUT_STEM_NAME + POPS_PREDS_SUFFIX),
            extension=".txt",
            read_spec=DataFrameReadSpec(
                DataFrameWhiteSpaceSepTextFormat(comment_code="#")
            ),
        )
        return cls(
            munge_task=munge_task,
            pops_run_task=pops_run_task,
            preds_result_table_task=preds_result_table_task,
        )
