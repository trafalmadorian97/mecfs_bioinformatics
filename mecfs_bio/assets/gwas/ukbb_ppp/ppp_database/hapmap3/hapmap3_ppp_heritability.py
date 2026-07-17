"""
Per-protein LDSC heritability over the HapMap3 UKB-PPP database: all-variants and
cis-excluded heritability for every protein, plus the per-protein sample-size table it
depends on.
"""

from mecfs_bio.asset_generator.ukbb_ppp_slim_protein_asset_generator import (
    generate_ppp_sample_size_task,
)
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.hapmap3_ppp_database_protein_files import (
    HAPMAP_3_PPP_DATABASE,
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

# Sample sizes are index-independent (a property of each protein's GWAS), so the asset id
# carries no index label.
PPP_SAMPLE_SIZES = generate_ppp_sample_size_task(asset_id="ppp_sample_size")

HAPMAP_3_PPP_HERITABILITY = PppProteinHeritabilityTask.create(
    asset_id="ppp_heritability_hapmap_3",
    protein_tasks=HAPMAP_3_PPP_DATABASE.protein_tasks,
    index_task=HAPMAP_3_PPP_DATABASE_INDEX,
    ld_ref_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
    sample_size_task=PPP_SAMPLE_SIZES,
    st3_task=SUN_ET_AL_2023_PQTL_SUPPLEMENTARY_RAW,
)
