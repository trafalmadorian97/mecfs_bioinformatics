"""
Sites-only allele-frequency tables of the 1000 Genomes EUR panels, for
genome-reference harmonization.

Wrapped in DiscardDepsWrapper so the multi-gigabyte genotype VCFs are not kept in the
asset store, only the derived parquet.
"""

from mecfs_bio.assets.reference_data.thousand_genomes.eur_hg19_vcf import (
    THOUSAND_GENOMES_EUR_HG19_VCF,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_hg38_30x_vcf import (
    THOUSAND_GENOMES_EUR_HG38_30X_VCF,
)
from mecfs_bio.build_system.task.discard_deps_task_wrapper import DiscardDepsWrapper
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    ReferencePanelAlleleFrequencyTask,
)

THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES = DiscardDepsWrapper(
    ReferencePanelAlleleFrequencyTask.create(
        vcf_task=THOUSAND_GENOMES_EUR_HG19_VCF,
        asset_id="thousand_genomes_eur_hg19_panel_allele_frequencies",
        build="19",
    )
)

THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES = DiscardDepsWrapper(
    ReferencePanelAlleleFrequencyTask.create(
        vcf_task=THOUSAND_GENOMES_EUR_HG38_30X_VCF,
        asset_id="thousand_genomes_eur_hg38_panel_allele_frequencies",
        build="38",
    )
)
