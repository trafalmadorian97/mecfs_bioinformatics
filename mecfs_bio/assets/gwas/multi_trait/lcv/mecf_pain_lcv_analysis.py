"""
Task to use the latent causal variable method to estimate the causal direction multisite pain and ME/CFS
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_standard_analysis import (
    JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genomes_phase_3_v1_consolidated import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
)
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_chrom_pos_alleles import (
    HarmonizeGWASWithReferenceViaAlleles,
)
from mecfs_bio.build_system.task.lcv.lcv_task import LCVConfig, LCVTask
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.filter_rows_by_min_in_col import (
    FilterRowsByMinInCol,
)
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
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

MECFS_HARMONIZE_WITH_PAIN = HarmonizeGWASWithReferenceViaAlleles.create(
    asset_id="decode_me_harmonized_with_johnston",
    gwas_data_task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS.join_task,
    reference_task=JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
    palindrome_strategy="drop",
    gwas_pipe=CompositePipe(
        [
            RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
            UniquePipe(
                by=[
                    GWASLAB_CHROM_COL,
                    GWASLAB_POS_COL,
                    GWASLAB_EFFECT_ALLELE_COL,
                    GWASLAB_NON_EFFECT_ALLELE_COL,
                ],
                keep="first",
                order_by=[GWASLAB_RSID_COL],
            ),
        ]
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


MECFS_PAIN_LCV_ANALYSIS = LCVTask.create(
    asset_id="mecfs_pain_lcv_analysis",
    trait_1_data=MECFS_HARMONIZE_WITH_PAIN,
    trait_2_data=JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
    consolidated_ld_scores=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
    config=LCVConfig(50),
    trait_1_pipe=CompositePipe([FilterRowsByMinInCol(1e-15, col=GWASLAB_SE_COL)]),
    trait_2_pipe=CompositePipe(
        [
            FilterRowsByMinInCol(1e-15, col=GWASLAB_SE_COL),
            ToPolarsPipe(),
        ]
    ),
)
