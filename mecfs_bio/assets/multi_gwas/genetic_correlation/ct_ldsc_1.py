from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_sumstats_liftover_to_37 import \
    DECODE_ME_GWAS_1_SUMSTATS_LIFTOVER_TO_37
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_with_annovar_37_rsids_sumstats import \
    DECODEME_ME_SUMSTATS_37_WITH_ANNOVAR_RSID
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_standard_analysis import \
    JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genome_phase_3_v1_extracted import \
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import GeneticCorrelationByCTLDSCTask, \
    SumstatsSource

CT_LDSC_INITIAL =GeneticCorrelationByCTLDSCTask.create(
    "initial_genetic_correlation_by_ct_ldsc",
    sources=[
SumstatsSource(
DECODEME_ME_SUMSTATS_37_WITH_ANNOVAR_RSID,
    alias="me_cfs"
),
        SumstatsSource(
    JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="multsite_pain"
        ),


    ],
    ld_ref_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED
)
