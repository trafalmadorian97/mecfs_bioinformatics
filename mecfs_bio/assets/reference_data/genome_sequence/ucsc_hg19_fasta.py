"""
UCSC hg19 (GRCh37) genome sequence, pinned by the md5 UCSC publishes in md5sum.txt.

The indexed form is what genome-reference harmonization reads. It is wrapped in
DiscardDepsWrapper, so only the uncompressed FASTA and its .fai are stored, not the
gzipped download as well. Note that UCSC hg19 chrM is not the rCRS mitochondrial
sequence, which is why harmonization excludes MT by default.
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

UCSC_HG19_FASTA_GZ = DownloadFileTask(
    meta=ReferenceFileMeta(
        id="ucsc_hg19_fasta_gz",
        group="genome_sequence",
        sub_group="ucsc_hg19",
        sub_folder=PurePath("raw"),
        extension=".fa.gz",
    ),
    url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
    md5_hash="806c02398f5ac5da8ffd6da2d1d5d1a9",
)

UCSC_HG19_INDEXED_FASTA = DiscardDepsWrapper(
    IndexedFastaTask.create(
        fasta_gz_task=UCSC_HG19_FASTA_GZ, asset_id="ucsc_hg19_indexed_fasta", build="19"
    )
)
