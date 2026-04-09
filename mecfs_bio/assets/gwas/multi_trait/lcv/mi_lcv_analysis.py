"""
Apply LCV to a list of a candidate upstream causal traits against the downstream trait Myocardial Infarction.

The goal here is mainly to make sure my LCV implementation is correct by testing for well known causal factors like LDL.
"""

from mecfs_bio.asset_generator.lcv_asset_generator import (
    LCVSourceTraitInfo,
    lcv_generate,
)
from mecfs_bio.assets.gwas.c_reactive_protein.said_et_al.analysis.said_crp_standard_analysis import \
    SAID_ET_AL_EUR_CRP_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.educational_attainment.lee_et_al_2018.processed_gwas_data.lee_et_al_magma_task_generator import \
    LEE_ET_AL_2018_COMBINED_MAGMA_TASKS
from mecfs_bio.assets.gwas.ldl.willer_et_al.analysis.willer_ldl_standard_analysis import (
    WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.myocardial_infarction.analysis.mi_standard_analysis import (
    MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.triglycerides.willer_et_al.analysis.triglycide_standard_analysis import (
    WILLER_ET_AL_EUR_TG_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genomes_phase_3_v1_consolidated import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
)
from mecfs_bio.build_system.task.lcv.lcv_task import LCVConfig
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaIfNeededPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe
from mecfs_bio.build_system.task.pipes.filter_rows_by_min_in_col import (
    FilterRowsByMinInCol,
)
from mecfs_bio.build_system.task.pipes.to_polars_pipe import ToPolarsPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_SE_COL

MI_LCV_TASK_GROUP = lcv_generate(
    base_name="mi_analysis",
    upstream_traits=[
        LCVSourceTraitInfo(
            name="LDL",
            df_task=WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
            pipe=CompositePipe([FilterRowsByMinInCol(1e-15, col=GWASLAB_SE_COL)]),
        ),
        LCVSourceTraitInfo(
            name="Triglycerides",
            df_task=WILLER_ET_AL_EUR_TG_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
            pipe=CompositePipe([FilterRowsByMinInCol(1e-15, col=GWASLAB_SE_COL)]),
        ),
        LCVSourceTraitInfo(
            name="CRP",
            df_task=SAID_ET_AL_EUR_CRP_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
            pipe=CompositePipe([FilterRowsByMinInCol(1e-15, col=GWASLAB_SE_COL)]),
        ),
        LCVSourceTraitInfo(
            name="Educational_Attainment",
            df_task=LEE_ET_AL_2018_COMBINED_MAGMA_TASKS.parquet_file_task,
            pipe=CompositePipe([FilterRowsByMinInCol(1e-15, col=GWASLAB_SE_COL)]),
        )
    ],
    downstream_traits=[
        LCVSourceTraitInfo(
            name="MI",
            df_task=MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
            pipe=CompositePipe(
                [
                    ComputeBetaIfNeededPipe(),
                    ComputeSEPipe(),
                    FilterRowsByMinInCol(1e-15, col=GWASLAB_SE_COL),
                    ToPolarsPipe(),
                ]
            ),
        )
    ],
    consolidated_ld_scores_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
    config=LCVConfig(chisq_exclude_factor_threshold=50),
)
