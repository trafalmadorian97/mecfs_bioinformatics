"""
HapMap3 regression-SNP index for fast all-pairs CT-LDSC: the GenomicSEM w_hm3 SNP list joined to
the 1000 Genomes phase 3 EUR LD scores, with hg19 positions.
"""

from mecfs_bio.assets.reference_data.genomic_sem_reference.hapmap3_snpist import (
    HAPMAP3_SNPLIST_FOR_GENOMIC_SEM,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genomes_phase_3_v1_consolidated import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
)
from mecfs_bio.build_system.task.fast_ct_ldsc.variant_index_task import (
    FastCtLdscVariantIndexTask,
)

HAPMAP3_EUR_LDSC_VARIANT_INDEX = FastCtLdscVariantIndexTask.create(
    "hapmap3_eur_ldsc_variant_index",
    hapmap3_task=HAPMAP3_SNPLIST_FOR_GENOMIC_SEM,
    ld_scores_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
)
