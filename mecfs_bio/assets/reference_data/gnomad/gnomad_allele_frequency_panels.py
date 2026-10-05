"""
gnomAD genome allele-frequency panels for genome-reference harmonization.

The hg19 panel is built from v2.1.1 (about 14 h streaming, about 3.7 GiB). The hg38 panel
from v4.1 is defined but built only on request (about 14 h, about 10 GiB). Neither is wrapped
in DiscardDepsWrapper; the per-chromosome parts stay in the asset store under
reference_data/gnomad/<release>/per_chromosome/, which a path_remap rule may move.
"""

from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.gnomad.gnomad_releases import (
    GNOMAD_V2_1_1_GENOMES,
    GNOMAD_V4_1_GENOMES,
)
from mecfs_bio.build_system.task_generator.gnomad_allele_frequency_panel_task_generator import (
    generate_gnomad_allele_frequency_panel_tasks,
)

GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES = (
    generate_gnomad_allele_frequency_panel_tasks(
        asset_id="gnomad_v2_1_1_genomes_hg19_allele_frequencies",
        release=GNOMAD_V2_1_1_GENOMES,
        fasta_task=UCSC_HG19_INDEXED_FASTA,
    ).panel_task
)

GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES = (
    generate_gnomad_allele_frequency_panel_tasks(
        asset_id="gnomad_v4_1_genomes_hg38_allele_frequencies",
        release=GNOMAD_V4_1_GENOMES,
        fasta_task=UCSC_HG38_INDEXED_FASTA,
    ).panel_task
)
