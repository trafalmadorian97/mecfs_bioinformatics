from typing import Sequence, Mapping

from attrs import frozen

from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.concat_frames_task import ConcatFramesTask
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_chrom_pos_alleles import \
    HarmonizeGWASWithReferenceViaAlleles
from mecfs_bio.build_system.task.lcv.lcv_task import LCVTask, LCVConfig
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.pipes.uniquepipe import UniquePipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL, \
    GWASLAB_NON_EFFECT_ALLELE_COL, GWASLAB_RSID_COL


@frozen
class LCVRun:
    harmonization_task:Task
    lcv_task:Task

@frozen
class LCVTaskGroup:
    lcv_run_mapping: Mapping[tuple[str,str],LCVRun]
    agg_task: Task


@frozen
class LCVSourceTraitInfo:
    name: str
    df_task: Task
    pipe: DataProcessingPipe = IdentityPipe()

def lcv_generate(
        base_name:str,
        upstream_traits: Sequence[LCVSourceTraitInfo],
        downstream_traits: Sequence[LCVSourceTraitInfo],
        consolidated_ld_scores_task: Task,
        config: LCVConfig
)-> LCVTaskGroup:
    unique_pipe=UniquePipe(
                            by=[
                                GWASLAB_CHROM_COL,
                                GWASLAB_POS_COL,
                                GWASLAB_EFFECT_ALLELE_COL,
                                GWASLAB_NON_EFFECT_ALLELE_COL,
                            ],
                            keep="first",
                            order_by=[GWASLAB_RSID_COL],
                        )
    run_mapping = {}
    task_list=[]
    for us in upstream_traits:
        for ds in downstream_traits:
            harmonization_task =HarmonizeGWASWithReferenceViaAlleles.create (
                asset_id=f"{base_name}_{us.name}_harmonize_with_{ds.name}_for_lcv",
                gwas_data_task=us.df_task,
                reference_task=ds.df_task,
                gwas_pipe=CompositePipe(
                    [
                        us.pipe,
                        unique_pipe,
                    ]
                ),
                ref_pipe=CompositePipe(
                    [
                        ds.pipe,
                        unique_pipe,
                    ]
                ),
                palindrome_strategy="drop"
            )

            lcv_task = LCVTask.create(
                 base_name+"_lcv_"+us.name+ "_ " + ds.name,
                    trait_1_data=harmonization_task,
                    trait_2_data=ds.df_task,
                    consolidated_ld_scores= consolidated_ld_scores_task,
                    config=config,
                    trait_2_pipe=CompositePipe(
                    [
                        ds.pipe,
                        unique_pipe,
                    ]
                )
                )
            run_mapping[(us.name, ds.name)] =LCVRun(
               lcv_task=lcv_task,
                harmonization_task=harmonization_task
            )
            task_list.append(lcv_task)
    agg_task = ConcatFramesTask.create(
        asset_id=base_name+"_lcv_agg",
        frames_tasks=task_list,
        out_format=ParquetOutFormat()
    )
    return LCVTaskGroup(
        lcv_run_mapping=run_mapping,
        agg_task=agg_task
    )



