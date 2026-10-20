"""
Task generator applying vectorized cell-type-specific S-LDSC to one GWAS across several prepared
annotation sets, with the same downstream tables and plots as SLDSCTaskGenerator. Agnostic to the
source of the annotation sets; every generated asset id carries the _vec suffix.
"""

from typing import Mapping, Sequence

from attrs import frozen

from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.multiple_testing_table_task import Method
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.sldsc.vectorized_cell_analysis_by_sldsc_task import (
    VectorizedCellAnalysisBySLDSCTask,
)
from mecfs_bio.build_system.task_generator.sldsc_task_generator import (
    VEC_ID_SUFFIX,
    CellAnalysisTaskGroup,
    CellOrTissueLabelRecord,
    build_cell_analysis_task_group,
)
from mecfs_bio.constants.sldsc_constants import (
    SLDSC_DEFAULT_MIN_HAPMAP3_SNPS,
    SLDSCWeighting,
)


@frozen(slots=True)
class PreparedAnnotationSetRecord:
    entry_name: str
    annotation_set_task: Task
    cell_or_tissue_labels_task: CellOrTissueLabelRecord | None


@frozen(slots=True)
class VectorizedSLDSCTaskGenerator:
    partitioned_tasks: Mapping[str, CellAnalysisTaskGroup]

    def get_terminal_tasks(self) -> Sequence[Task]:
        result: list[Task] = []
        for group in self.partitioned_tasks.values():
            result.extend(group.terminal_tasks())
        return result

    @classmethod
    def create(
        cls,
        base_name: str,
        source_sumstats_task: Task,
        baseline_task: Task,
        entries: Sequence[PreparedAnnotationSetRecord],
        multiple_testing_alpha: float,
        multiple_testing_method: Method,
        weighting: SLDSCWeighting = "exact",
        pre_pipe: DataProcessingPipe = IdentityPipe(),
        min_hapmap3_snps: int = SLDSC_DEFAULT_MIN_HAPMAP3_SNPS,
    ) -> "VectorizedSLDSCTaskGenerator":
        return cls(
            {
                entry.entry_name: build_cell_analysis_task_group(
                    base_name=base_name,
                    entry_name=entry.entry_name,
                    cell_task=VectorizedCellAnalysisBySLDSCTask.create(
                        asset_id=base_name
                        + "_"
                        + entry.entry_name
                        + "_cell_analysis_by_s_ldsc"
                        + VEC_ID_SUFFIX,
                        source_sumstats_task=source_sumstats_task,
                        baseline_task=baseline_task,
                        annotation_set_task=entry.annotation_set_task,
                        weighting=weighting,
                        pre_pipe=pre_pipe,
                        min_hapmap3_snps=min_hapmap3_snps,
                    ),
                    labels=entry.cell_or_tissue_labels_task,
                    multiple_testing_alpha=multiple_testing_alpha,
                    multiple_testing_method=multiple_testing_method,
                    id_suffix=VEC_ID_SUFFIX,
                )
                for entry in entries
            }
        )
