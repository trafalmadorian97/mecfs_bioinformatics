"""
Pan-UK Biobank allele-frequency panel (GRCh37, imputed UK Biobank, six genetic-ancestry
groups) for genome-reference harmonization.

295 manifest rows (chr21 47, chr22 86, X 162) have ref and alt swapped relative to the hg19
FASTA, inherited from UK Biobank's imputed BGEN allele order; the panel drops them and the
build fails if the count changes. Measured by experiments/claude/pan_ukbb_manifest/.
"""

from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.pan_ukbb.pan_ukbb_variant_manifest import (
    PAN_UKBB_VARIANT_MANIFEST,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.pan_ukbb.pan_ukbb_allele_frequency_panel_task import (
    PanUkbbAlleleFrequencyPanelTask,
)

PAN_UKBB_HG19_ALLELE_FREQUENCIES = PanUkbbAlleleFrequencyPanelTask.create(
    asset_id="pan_ukbb_hg19_allele_frequencies",
    manifest_task=PAN_UKBB_VARIANT_MANIFEST,
    fasta_task=UCSC_HG19_INDEXED_FASTA,
    groups=("ukb_afr", "ukb_amr", "ukb_csa", "ukb_eas", "ukb_eur", "ukb_mid"),
    chromosomes=(*range(1, 23), 23),
    expected_ref_mismatches=295,
)
