from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genome_phase_3_v1_extracted import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
)
from mecfs_bio.build_system.task.consolidate_ld_scores_task import (
    ConsolidateLDScoresTask,
)

THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE = ConsolidateLDScoresTask.create(
    asset_id="thousand_genomes_phase_3_v1_eur_ld_scores_consolidated",
    extracted_ld_score_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
)
