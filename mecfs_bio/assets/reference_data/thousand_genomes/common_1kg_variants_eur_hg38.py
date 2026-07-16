"""
Common biallelic SNVs from the 1000 Genomes 30x EUR panel (hg38), chrX included.

The membership list behind the common-1kg PPP variant index mode. CHR is coded the
way gwaslab codes it, so X appears as 23.
"""

from mecfs_bio.assets.reference_data.thousand_genomes.eur_hg38_30x_vcf import (
    THOUSAND_GENOMES_EUR_HG38_30X_VCF,
)
from mecfs_bio.build_system.task.ppp_database.filter_common_1kg_variants_task import (
    FilterCommon1kgVariantsTask,
)

COMMON_1KG_VARIANTS_EUR_HG38 = FilterCommon1kgVariantsTask.create(
    thousand_genomes=THOUSAND_GENOMES_EUR_HG38_30X_VCF,
    asset_id="common_1kg_variants_eur_hg38",
)
