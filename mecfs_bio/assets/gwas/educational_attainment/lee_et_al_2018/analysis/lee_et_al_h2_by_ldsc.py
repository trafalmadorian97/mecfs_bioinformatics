from mecfs_bio.assets.gwas.educational_attainment.lee_et_al_2018.processed_gwas_data.lee_et_al_magma_task_generator import (
    LEE_ET_AL_2018_COMBINED_MAGMA_TASKS,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genome_phase_3_v1_extracted import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    QuantPhenotype,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_snp_heritability_by_ldsc_task import (
    SNPHeritabilityByLDSCTask,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe

LEE_ET_AL_2018_H2_BY_LDSC = SNPHeritabilityByLDSCTask.create(
    asset_id="lee_et_al_2018_r2_by_ldsc",
    source_sumstats_task=LEE_ET_AL_2018_COMBINED_MAGMA_TASKS.sumstats_task,
    ld_ref_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
    set_sample_size=257841,  # source: 30038396-GCST006572-EFO_0008354.h.tsv.gz-meta.yaml
    pipe=IdentityPipe(),
    phenotype_info=QuantPhenotype(),
    build="19",
)
