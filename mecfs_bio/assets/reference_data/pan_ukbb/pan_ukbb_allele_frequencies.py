"""
Pan-UK Biobank allele-frequency panel (GRCh37, imputed UK Biobank, six genetic-ancestry
groups) for genome-reference harmonization.

On chr21, chr22 and X, genotyped manifest rows without a gnomAD frequency have ref and alt
swapped; the Task drops them (drop_swap_profile_rows). They include all 295 rows whose ref
disagrees with the hg19 FASTA (chr21 47, chr22 86, X 162; measured by
experiments/claude/pan_ukbb_manifest/), so no FASTA mismatch remains and the build fails if
one appears.
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
    expected_ref_mismatches=0,
)
