"""
Task to download UKBB LD matrix for a small subset of the genome
Chromosome 1 from 174000001 to 17700000 (genome build GRCh37).
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

CHR1_173000001_17600000_UKBB_LD_MATRIX_DOWNLOAD = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="ukbb_reference_ld",
        sub_group="chr1_173000001_17600000",
        sub_folder=PurePath("raw"),
        id=AssetId("ukbb_chr1_173000001_17600000_ld_matrix"),
        filename="chr1_173000001_176000001",
        extension=".npz",
    ),
    url="https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/chr1_173000001_176000001.npz",
    md5_hash=None,
)
