"""
UCSC hg38 (GRCh38) genome sequence, pinned by the md5 UCSC publishes in md5sum.txt.

Used to tune the stringent ambiguous-indel rules on build-38 DecodeME, whose source
orientation is fully reference-consistent. Wrapped in DiscardDepsWrapper so only the
uncompressed FASTA and its .fai are stored.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.discard_deps_task_wrapper import DiscardDepsWrapper
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
from mecfs_bio.build_system.task.genome_reference_harmonization.indexed_fasta_task import (
    IndexedFastaTask,
)

UCSC_HG38_FASTA_GZ = DownloadFileTask(
    meta=ReferenceFileMeta(
        id="ucsc_hg38_fasta_gz",
        group="genome_sequence",
        sub_group="ucsc_hg38",
        sub_folder=PurePath("raw"),
        extension=".fa.gz",
    ),
    url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
    md5_hash="1c9dcaddfa41027f17cd8f7a82c7293b",
)

UCSC_HG38_INDEXED_FASTA = DiscardDepsWrapper(
    IndexedFastaTask.create(
        fasta_gz_task=UCSC_HG38_FASTA_GZ, asset_id="ucsc_hg38_indexed_fasta"
    )
)
