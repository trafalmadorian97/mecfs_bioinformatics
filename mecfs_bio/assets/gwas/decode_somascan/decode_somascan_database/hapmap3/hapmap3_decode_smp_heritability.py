"""
Per-protein LDSC heritability over the HapMap3 deCODE SomaScan SMP database: all-variants and
cis-excluded heritability for every assay, with the 1000 Genomes European LD-score reference
(the Icelandic cohort is of European ancestry) and gene coordinates from the MAGMA GRCh38 table.
"""

from mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_database_protein_files import (
    HAPMAP_3_DECODE_SMP_DATABASE,
)
from mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_index import (
    HAPMAP_3_DECODE_SMP_INDEX,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genomes_phase_3_v1_consolidated import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.processed.magma_entrez_gene_locations_data_build_38_unzipped import (
    MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_38_EXTRACTED,
)
from mecfs_bio.build_system.task.ppp_ldsc.decode_protein_heritability_task import (
    DecodeProteinHeritabilityTask,
)

HAPMAP_3_DECODE_SMP_HERITABILITY = DecodeProteinHeritabilityTask.create(
    asset_id="decode_somascan_smp_heritability_hapmap_3",
    protein_tasks=HAPMAP_3_DECODE_SMP_DATABASE.protein_tasks,
    index_task=HAPMAP_3_DECODE_SMP_INDEX,
    ld_scores_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
    gene_coords_task=MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_38_EXTRACTED,
)
