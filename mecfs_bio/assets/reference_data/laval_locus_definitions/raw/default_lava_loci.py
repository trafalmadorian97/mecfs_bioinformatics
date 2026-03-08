from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

LAVA_DEFAULT_LOCUS_REFERENCE = DownloadFileTask(
    url="https://raw.githubusercontent.com/josefin-werme/LAVA/refs/heads/main/support_data/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile",
    md5_hash=None,
    meta=ReferenceFileMeta(
        group="lava",
        sub_group="default",
        sub_folder=PurePath("rar"),
        id="lava_default_locus_reference",
    ),
)
