"""
gnomAD genome allele-frequency panels for genome-reference harmonization.

The hg19 panel is built from v2.1.1 (about 14 h streaming, about 3.9 GiB). The hg38 panel
from v4.1 (about 14 h, about 10 GiB).

GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES_REHOSTED downloads a copy of the built hg19
panel instead, so that CI and new machines need not run the 14 h build. Its md5 is that of
the panel built on 2026-10-06, and its metadata (build, ancestry columns) is the built
panel's.
"""

from pathlib import PurePath

import attrs

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
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
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

GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES_REHOSTED = DownloadFileTask(
    meta=attrs.evolve(
        GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES.meta,
        id="gnomad_v2_1_1_genomes_hg19_allele_frequencies_rehosted",
        sub_folder=PurePath("rehosted"),
    ),
    url=(
        "https://www.dropbox.com/scl/fi/ccvdhj47qj454qcdndf3h/"
        "gnomad_allele_frequencies.parquet?rlkey=8qqhj9unoskcgmwcjfd4x40ee&dl=1"
    ),
    md5_hash="f10c58ea238ec758ac4c6754e8b73cce",
)

GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES = (
    generate_gnomad_allele_frequency_panel_tasks(
        asset_id="gnomad_v4_1_genomes_hg38_allele_frequencies",
        release=GNOMAD_V4_1_GENOMES,
        fasta_task=UCSC_HG38_INDEXED_FASTA,
    ).panel_task
)
