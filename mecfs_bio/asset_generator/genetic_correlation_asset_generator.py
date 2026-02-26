import itertools
from typing import Sequence, Mapping
from frozendict import frozendict
from attrs import frozen

from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genome_phase_3_v1_extracted import \
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.concat_frames_task import ConcatFramesTask
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import GenomeBuild
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import SumstatsSource, \
    GeneticCorrelationByCTLDSCTask
from mecfs_bio.build_system.task.pipe_dataframe_task import CSVOutFormat


@frozen
class GeneticCorrTasks:
    corr_tasks: Mapping[str, Task]
    aggregation_task: Task

def genetic_corr_by_ct_ldsc_asset_generator(
        base_name: str,
        sources: Sequence[SumstatsSource],
        ld_ref_task:Task= THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
        ld_file_filename_pattern: str = "/LDscore.@",
        build:GenomeBuild="19",

) ->GeneticCorrTasks:
    tasks={}
    for source_1, source_2 in itertools.combinations(sources, 2):
        sorted_sources = sorted([source_1, source_2], key=lambda x: x.alias.lower())
        task_name  = (base_name + "_ct_ldsc_corr_"+sorted_sources[0].alias+"__"+sorted_sources[1].alias).lower()
        task= GeneticCorrelationByCTLDSCTask.create(
            asset_id=task_name,
            sources=sorted_sources,
            ld_ref_task=ld_ref_task,
            build=build,
             ld_file_filename_pattern=ld_file_filename_pattern,
        )
        tasks[task_name] = task
    aggregation = ConcatFramesTask.create(
        asset_id=base_name+"_ct_ldsc_corr_aggregation",
        frames_tasks=list(tasks.values()),
        out_format=CSVOutFormat(",")
    )

    return GeneticCorrTasks(
        corr_tasks=frozendict(tasks),
        aggregation_task=aggregation,
    )
