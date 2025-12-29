from attrs import frozen

from mecfs_bio.asset_generator.concrete_magma_asset_generator import concrete_magma_assets_generate
from mecfs_bio.asset_generator.concrete_sldsc_generator import standard_sldsc_task_generator
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import ValidGwaslabFormat
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task_generator.magma_task_generator import MagmaTaskGeneratorFromRaw
from mecfs_bio.build_system.task_generator.sldsc_task_generator import SLDSCTaskGenerator


@frozen
class StandardAnalysisTaskGroup:
    sldsc_tasks:SLDSCTaskGenerator
    magma_tasks: MagmaTaskGeneratorFromRaw


def concrete_standard_analysis_generator(
        base_name: str,
        raw_gwas_data_task: Task,
        fmt: ValidGwaslabFormat,
        sample_size: int,
        pre_pipe: DataProcessingPipe = IdentityPipe(),
        pre_sldsc_pipe: DataProcessingPipe = IdentityPipe()
    )-> StandardAnalysisTaskGroup:
   magma_tasks =  concrete_magma_assets_generate(
       base_name=base_name,
        raw_gwas_data_task=raw_gwas_data_task,
       fmt=fmt,
       sample_size=sample_size,
       pre_pipe=pre_pipe,
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