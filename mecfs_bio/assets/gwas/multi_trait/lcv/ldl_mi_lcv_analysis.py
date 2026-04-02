from mecfs_bio.assets.gwas.ldl.million_veterans.analysis.million_veterans_ldl_eur_magma_task_generator import (
    MILLION_VETERANS_EUR_LDL_MAGMA_TASKS,
)
from mecfs_bio.assets.gwas.myocardial_infarction.analysis.mi_standard_analysis import (
    MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genomes_phase_3_v1_consolidated import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
)
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_chrom_pos_alleles import (
    HarmonizeGWASWithReferenceViaAlleles,
)
from mecfs_bio.build_system.task.lcv.lcv_task import LCVConfig, LCVTask
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaIfNeededPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe
from mecfs_bio.build_system.task.pipes.filter_rows_by_min_in_col import (
    FilterRowsByMinInCol,
)
from mecfs_bio.build_system.task.pipes.to_polars_pipe import ToPolarsPipe
from mecfs_bio.build_system.task.pipes.uniquepipe import UniquePipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SE_COL,
)

LDL_HARMONIZE_WITH_MI = HarmonizeGWASWithReferenceViaAlleles.create(
    asset_id="ldl_harmonized_with_mi",
    gwas_data_task=MILLION_VETERANS_EUR_LDL_MAGMA_TASKS.parquet_file_task,
    reference_task=MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
    palindrome_strategy="drop",
    gwas_pipe=UniquePipe(
        by=[
            GWASLAB_CHROM_COL,
            GWASLAB_POS_COL,
            GWASLAB_EFFECT_ALLELE_COL,
            GWASLAB_NON_EFFECT_ALLELE_COL,
        ],
        keep="first",
        order_by=[GWASLAB_RSID_COL],
    ),
    ref_pipe=UniquePipe(
        by=[
            GWASLAB_CHROM_COL,
            GWASLAB_POS_COL,
            GWASLAB_EFFECT_ALLELE_COL,
            GWASLAB_NON_EFFECT_ALLELE_COL,
        ],
        keep="first",
        order_by=[GWASLAB_RSID_COL],
    ),
)


LDL_MI_LCV_ANALYSIS = LCVTask.create(
    asset_id="ldl_mi_lcv_analysis",
    trait_1_data=LDL_HARMONIZE_WITH_MI,
    trait_2_data=MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
    consolidated_ld_scores=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
    config=LCVConfig(50),
    trait_1_pipe=CompositePipe([FilterRowsByMinInCol(1e-15, col=GWASLAB_SE_COL)]),
    trait_2_pipe=CompositePipe(
        [
            ComputeBetaIfNeededPipe(),
            ComputeSEPipe(),
            FilterRowsByMinInCol(1e-15, col=GWASLAB_SE_COL),
            ToPolarsPipe(),
        ]
    ),
)
