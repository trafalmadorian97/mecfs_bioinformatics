"""
Download the default LAVA partitioning of the genome into loci.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

DEFAULT_LAVA_LOCUS_FILE = DownloadFileTask(
    ReferenceFileMeta(
        group="lava_locus_file",
        sub_group="default",
        sub_folder=PurePath("raw"),
        id="lava_default_locus_file",
        filename="blocks_s2500_m25_f1_w200.GRCh37_hg19",
        extension=".locfile",
    ),
    url="https://raw.githubusercontent.com/josefin-werme/LAVA/refs/heads/main/support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile",
    md5_hash=None,
)
