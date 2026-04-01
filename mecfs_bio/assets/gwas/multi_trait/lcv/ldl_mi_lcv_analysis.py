from mecfs_bio.assets.gwas.ldl.million_veterans.analysis.million_veterans_ldl_eur_magma_task_generator import (
    MILLION_VETERANS_EUR_LDL_MAGMA_TASKS,
)
from mecfs_bio.assets.gwas.myocardial_infarction.analysis.mi_standard_analysis import (
    MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genomes_phase_3_v1_consolidated import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
)
from mecfs_bio.build_system.task.lcv.lcv_task import LCVConfig, LCVTask
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaIfNeededPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe

LDL_MI_LCV_ANALYSIS = LCVTask.create(
    asset_id="ldl_mi_lcv_analysis",
    trait_1_data=MILLION_VETERANS_EUR_LDL_MAGMA_TASKS.parquet_file_task,
    trait_2_data=MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
    consolidated_ld_scores=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
    config=LCVConfig(50),
    trait_2_pipe=CompositePipe(
        [
ComputeBetaIfNeededPipe(),
            ComputeSEPipe()
        ]
    )
)
