from attrs import frozen

from mecfs_bio.asset_generator.annovar_37_basic_rsid_assignment import annovar_37_basic_rsid_assignment, \
    RSIDAssignmentTaskGroup
from mecfs_bio.asset_generator.concrete_magma_asset_generator import concrete_magma_assets_generate
from mecfs_bio.asset_generator.concrete_sldsc_generator import standard_sldsc_task_generator
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_constants import GWASLAB_SAMPLE_SIZE_COLUMN
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import ValidGwaslabFormat, \
    GWASLabCreateSumstatsTask
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe
from mecfs_bio.build_system.task_generator.magma_task_generator import MagmaTaskGeneratorFromRaw
from mecfs_bio.build_system.task_generator.sldsc_task_generator import SLDSCTaskGenerator


@frozen
class StandardAnalysisTaskGroup:
    sldsc_tasks:SLDSCTaskGenerator
    magma_tasks: MagmaTaskGeneratorFromRaw
    def get_terminal_tasks(self) ->list[Task]:
        return list(self.sldsc_tasks.get_terminal_tasks()) + self.magma_tasks.inner.terminal_tasks()


def concrete_standard_analysis_generator_assumme_already_has_rsid(
        base_name: str,
        raw_gwas_data_task: Task,
        fmt: ValidGwaslabFormat,
        sample_size: int,
        sample_size_for_sldsc: int|None = None,
        pre_pipe: DataProcessingPipe = IdentityPipe(),
        pre_sldsc_pipe: DataProcessingPipe = IdentityPipe()
    )-> StandardAnalysisTaskGroup:
   """
   Generate standard MAGMA and S-LDSC analysis tasks for given GWAS data,
   assuming that gwas data already contains rsids
   """
   magma_tasks =  concrete_magma_assets_generate(
       base_name=base_name,
        raw_gwas_data_task=raw_gwas_data_task,
       fmt=fmt,
       sample_size=sample_size,
       pre_pipe=pre_pipe,
   )
   sldsc_ss = sample_size_for_sldsc if sample_size_for_sldsc is not None else sample_size
   pre_sldsc_pipe= CompositePipe(
       [
           SetColToConstantPipe(col_name=GWASLAB_SAMPLE_SIZE_COLUMN, constant=sldsc_ss),
           pre_sldsc_pipe
       ]
   )

   sldsc_tasks= standard_sldsc_task_generator(
       sumstats_task=magma_tasks.sumstats_task,
       base_name=base_name,
       pre_pipe=pre_sldsc_pipe,
   )
   return StandardAnalysisTaskGroup(
       sldsc_tasks=sldsc_tasks,
       magma_tasks=magma_tasks,
   )



@frozen
class StandardAnalysisTaskGroupAddRSIDS:
    tasks: StandardAnalysisTaskGroup
    initial_sumstats_task: Task
    assign_rsids_task_group: RSIDAssignmentTaskGroup
    def terminal_tasks(self) -> list[Task]:
        return list(self.tasks.get_terminal_tasks())

def concrete_standard_analysis_generator_no_rsid(
        base_name: str,
        raw_gwas_data_task: Task,
        fmt: ValidGwaslabFormat,
        sample_size: int,
        sample_size_for_sldsc: int|None = None,
        pre_pipe: DataProcessingPipe = IdentityPipe(),
        pre_sldsc_pipe: DataProcessingPipe = IdentityPipe()
)-> StandardAnalysisTaskGroupAddRSIDS:
    sumstats_37_task = GWASLabCreateSumstatsTask(
        df_source_task=raw_gwas_data_task,
        asset_id=AssetId(base_name + "_initial_sumstats_37"),
        basic_check=True,
        genome_build="infer",
        liftover_to="19",
        fmt=fmt,
        pre_pipe=pre_pipe,
    )
    rsids_assigned_task_group = annovar_37_basic_rsid_assignment(
        sumstats_task=sumstats_37_task,
        base_name=base_name,
        use_gwaslab_rsids_convention=True
    )
    standard_tasks = concrete_standard_analysis_generator_assumme_already_has_rsid(
        base_name=base_name,
        raw_gwas_data_task=rsids_assigned_task_group.join_task,
        fmt="gwaslab",
        sample_size=sample_size,
        pre_pipe=pre_pipe,
        pre_sldsc_pipe=pre_sldsc_pipe,
        sample_size_for_sldsc=sample_size_for_sldsc,
    )
    return StandardAnalysisTaskGroupAddRSIDS(
        tasks=standard_tasks,
        initial_sumstats_task=sumstats_37_task,
        assign_rsids_task_group=rsids_assigned_task_group,
    )

