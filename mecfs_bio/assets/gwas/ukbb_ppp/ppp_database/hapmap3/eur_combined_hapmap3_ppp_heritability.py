"""
Per-protein LDSC heritability over the HapMap3 UKB-PPP database for the full European
(combined) cohort (discovery + replication, n = 52,363): all-variants and cis-excluded
heritability for every protein.

Parallel to eur_discovery_hapmap3_ppp_heritability, but the combined slim protein files
store the per-protein sample size N (include_sample_size defaults on), so no separate
ranged-read sample-size table is needed: sample_size_task=None makes the heritability task
read the constant N straight from each slim parquet. Distinct asset ids let it coexist with
the discovery assets.
"""

from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.eur_combined_hapmap3_ppp_database_protein_files import (
    HAPMAP_3_PPP_DATABASE_EUR_COMBINED,
)
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.hapmap_3_ppp_index import (
    HAPMAP_3_PPP_DATABASE_INDEX,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genome_phase_3_v1_extracted import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
)
from mecfs_bio.assets.reference_data.pqtls.raw.sun_et_al_2023_pqtls import (
    SUN_ET_AL_2023_PQTL_SUPPLEMENTARY_RAW,
)
from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_heritability_task import (
    PppProteinHeritabilityTask,
)

HAPMAP_3_PPP_HERITABILITY_EUR_COMBINED = PppProteinHeritabilityTask.create(
    asset_id="ppp_heritability_hapmap_3_eur_combined",
    protein_tasks=HAPMAP_3_PPP_DATABASE_EUR_COMBINED.protein_tasks,
    index_task=HAPMAP_3_PPP_DATABASE_INDEX,
    ld_ref_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
    sample_size_task=None,
    st3_task=SUN_ET_AL_2023_PQTL_SUPPLEMENTARY_RAW,
)
