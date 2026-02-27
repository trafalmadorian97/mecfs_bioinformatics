from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

CHR1_174000001_17700000_UKBB_LD_LABELS_DOWNLOAD = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="ukbb_reference_ld",
        sub_group="chr1_174000001_17700000",
        sub_folder=PurePath("raw"),
        id=AssetId("ukbb_chr1_174000001_17700000_ld_labels"),
        filename="chr1_174000001_177000001",
        extension=".gz",
    ),
    url="https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/chr1_174000001_177000001.gz",
    md5_hash=None,
)
