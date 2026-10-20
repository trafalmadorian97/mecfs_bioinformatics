"""
Asset generator applying vectorized S-LDSC with the six standard LDSC-SEG annotation sets
(Finucane et al. 2018), mirroring standard_sldsc_task_generator.
"""

from mecfs_bio.asset_generator.concrete_sldsc_generator import (
    IMMGEN_LABELS,
    MULTI_TISSUE_CHROMATIN_LABELS,
    MULTI_TISSUE_GENE_EXPR_LABELS,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.prepared.prepared_sldsc_reference_sets import (
    PREPARED_BASELINE_V1_2_HM3_WEIGHTS,
    PREPARED_CAHOY_ANNOTATION_SET,
    PREPARED_CORCES_ATAC_ANNOTATION_SET,
    PREPARED_GTEX_BRAIN_ANNOTATION_SET,
    PREPARED_IMMGEN_ANNOTATION_SET,
    PREPARED_MULTI_TISSUE_CHROMATIN_ANNOTATION_SET,
    PREPARED_MULTI_TISSUE_GENE_EXPR_ANNOTATION_SET,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task_generator.vectorized_sldsc_task_generator import (
    PreparedAnnotationSetRecord,
    VectorizedSLDSCTaskGenerator,
)
from mecfs_bio.constants.sldsc_constants import SLDSCWeighting


def standard_vectorized_sldsc_task_generator(
    sumstats_task: Task,
    base_name: str,
    pre_pipe: DataProcessingPipe = IdentityPipe(),
    weighting: SLDSCWeighting = "exact",
) -> VectorizedSLDSCTaskGenerator:
    return VectorizedSLDSCTaskGenerator.create(
        base_name=base_name,
        source_sumstats_task=sumstats_task,
        baseline_task=PREPARED_BASELINE_V1_2_HM3_WEIGHTS,
        entries=[
            PreparedAnnotationSetRecord(
                entry_name="cahoy_cns",
                annotation_set_task=PREPARED_CAHOY_ANNOTATION_SET,
                cell_or_tissue_labels_task=None,
            ),
            PreparedAnnotationSetRecord(
                entry_name="gtex_brain",
                annotation_set_task=PREPARED_GTEX_BRAIN_ANNOTATION_SET,
                cell_or_tissue_labels_task=None,
            ),
            PreparedAnnotationSetRecord(
                entry_name="immgen",
                annotation_set_task=PREPARED_IMMGEN_ANNOTATION_SET,
                cell_or_tissue_labels_task=IMMGEN_LABELS,
            ),
            PreparedAnnotationSetRecord(
                entry_name="corces_atac",
                annotation_set_task=PREPARED_CORCES_ATAC_ANNOTATION_SET,
                cell_or_tissue_labels_task=None,
            ),
            PreparedAnnotationSetRecord(
                entry_name="multi_tissue_chromatin",
                annotation_set_task=PREPARED_MULTI_TISSUE_CHROMATIN_ANNOTATION_SET,
                cell_or_tissue_labels_task=MULTI_TISSUE_CHROMATIN_LABELS,
            ),
            PreparedAnnotationSetRecord(
                entry_name="multi_tissue_gene_expression",
                annotation_set_task=PREPARED_MULTI_TISSUE_GENE_EXPR_ANNOTATION_SET,
                cell_or_tissue_labels_task=MULTI_TISSUE_GENE_EXPR_LABELS,
            ),
        ],
        multiple_testing_alpha=0.01,
        multiple_testing_method="fdr_bh",
        weighting=weighting,
        pre_pipe=pre_pipe,
    )
