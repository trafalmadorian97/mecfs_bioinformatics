import itertools
from typing import Sequence, Mapping

from attrs import frozen

from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genome_phase_3_v1_extracted import \
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import SumstatsSource, \
    GeneticCorrelationByCTLDSCTask


@frozen
class GeneticCorrTasks:
    corr_tasks: Mapping[frozenset, Task]
    aggregation_task: Task

def genetic_corr_by_ct_ldsc_asset_generator(
        base_name: str,
        sources: Sequence[SumstatsSource],
        ld_ref_task:Task= THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,

        ld_file_filename_pattern: str = "/LDscore.@",

) ->GeneticCorrTasks:
    tasks={}
    for source_1, source_2 in itertools.combinations(sources, 2):
        sorted_sources = sorted([source_1, source_2], key=lambda x: x.alias)
        task_name  = base_name + "_ct_ldsc_corr_"+sorted_sources[0].alias+"_"+sorted_sources[1].alias
        task= GeneticCorrelationByCTLDSCTask.create(
            asset_id=task_name,
            sources=sorted_sources,
            ld_ref_task=ld_ref_task,
            build="19",
             ld_file_filename_pattern=ld_file_filename_pattern,
        )
        tasks[frozenset(sorted_sources)] = task
