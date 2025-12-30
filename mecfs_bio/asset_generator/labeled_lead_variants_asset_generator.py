from attrs import frozen

from mecfs_bio.assets.reference_data.ensembl_biomart.gene_thesaurus import (
    GENE_THESAURUS,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
    ValidGwaslabFormat,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_lead_variants_task import (
    GwasLabLeadVariantsTask,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe


@frozen
class LabelLeadVariantsTasks:
    raw_sumstats_task: Task
    lead_variants_task: Task
    labeled_lead_variants_task: Task


def generate_tasks_labeled_lead_variants(
    base_name: str,
    raw_gwas_data_task: Task,
    fmt: ValidGwaslabFormat,
    pre_pipe: DataProcessingPipe = IdentityPipe(),
) -> LabelLeadVariantsTasks:
    raw_sumstats_task = GWASLabCreateSumstatsTask(
        df_source_task=raw_gwas_data_task,
        fmt=fmt,
        basic_check=True,
        genome_build="infer",
        pre_pipe=pre_pipe,
        asset_id=AssetId(base_name + "_raw_sumstats_task_no_liftover"),
    )
    lead_variants_task = GwasLabLeadVariantsTask(
        sumstats_task=raw_sumstats_task, short_id=base_name + "_gwaslab_lead_variants"
    )
    labeled_lead_variants_task = JoinDataFramesTask.create_from_result_df(
        asset_id=base_name + "_labeled_lead_variants",
        result_df_task=lead_variants_task,
        reference_df_task=GENE_THESAURUS,
        how="left",
        left_on=["GENE"],
        right_on=["Gene name"],
        df_1_pipe=SelectColPipe(["GENE"]),
        df_2_pipe=SelectColPipe(
            ["Gene name", "NCBI gene (formerly Entrezgene) ID", "Gene stable ID"]
        ),
    )
    return LabelLeadVariantsTasks(
        raw_sumstats_task=raw_sumstats_task,
        lead_variants_task=lead_variants_task,
        labeled_lead_variants_task=labeled_lead_variants_task,
    )
