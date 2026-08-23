"""
LDSC heritability for each aptamer in the Western et al. 2024 CSF proteomics database.

"""

from mecfs_bio.assets.gwas.csf_pqtl.csf_database.hapmap3.hapmap3_csf_database_aptamer_files import (
    HAPMAP_3_CSF_DATABASE,
)
from mecfs_bio.assets.gwas.csf_pqtl.csf_database.hapmap3.hapmap3_csf_index import (
    HAPMAP_3_CSF_DATABASE_INDEX,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genomes_phase_3_v1_consolidated import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
)
from mecfs_bio.build_system.task.csf_ldsc.csf_protein_heritability_task import (
    CsfProteinHeritabilityTask,
)

HAPMAP_3_CSF_HERITABILITY = CsfProteinHeritabilityTask.create(
    asset_id="csf_heritability_hapmap_3",
    aptamer_tasks=HAPMAP_3_CSF_DATABASE.aptamer_tasks,
    index_task=HAPMAP_3_CSF_DATABASE_INDEX,
    ld_scores_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
)
