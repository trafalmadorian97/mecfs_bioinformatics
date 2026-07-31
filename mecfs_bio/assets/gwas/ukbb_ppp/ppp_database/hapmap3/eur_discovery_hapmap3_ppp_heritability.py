"""
Per-protein LDSC heritability over the HapMap3 UKB-PPP database for the European
(discovery) cohort (n = 34,557): all-variants and cis-excluded heritability for every
protein, plus the per-protein sample-size table it depends on.

Pinned to the discovery cohort (manifest passed explicitly) to stay consistent with the
slim files already built for it. This is the European-ancestry cohort whose LD structure
matches the European HapMap3 LD-score reference LDSC uses.
"""

from mecfs_bio.asset_generator.ukbb_ppp_slim_protein_asset_generator import (
    EUR_DISCOVERY_PPP_MANIFEST_PATH,
    generate_ppp_sample_size_task,
)
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.eur_discovery_hapmap3_ppp_database_protein_files import (
    HAPMAP_3_PPP_DATABASE,
)
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.hapmap_3_ppp_index import (
    HAPMAP_3_PPP_DATABASE_INDEX,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genomes_phase_3_v1_consolidated import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
)
from mecfs_bio.assets.reference_data.pqtls.processed.sun_et_al_2023_st3_extracted import (
    SUN_ET_AL_2023_ST3_EXTRACTED,
)
from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_heritability_task import (
    PppProteinHeritabilityTask,
)

# Sample sizes are index-independent (a property of each protein's GWAS), so the asset id
# carries no index label.
PPP_SAMPLE_SIZES = generate_ppp_sample_size_task(
    asset_id="ppp_sample_size",
    manifest_path=EUR_DISCOVERY_PPP_MANIFEST_PATH,
)

HAPMAP_3_PPP_HERITABILITY = PppProteinHeritabilityTask.create(
    asset_id="ppp_heritability_hapmap_3",
    protein_tasks=HAPMAP_3_PPP_DATABASE.protein_tasks,
    index_task=HAPMAP_3_PPP_DATABASE_INDEX,
    ld_scores_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
    sample_size_task=PPP_SAMPLE_SIZES,
    gene_coords_task=SUN_ET_AL_2023_ST3_EXTRACTED,
)
