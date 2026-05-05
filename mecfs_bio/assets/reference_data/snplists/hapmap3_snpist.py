"""
Genomic SEM requires a reference file containing a list of SNPs to use for ldsc.
Usually, these consist of the hapmap3 SNPs

The Task here downloads the reference file recommended by the GenomicSEM tutorial.  See"
https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Genome%E2%80%90wide-Models

I have rehosted the file on dropbox.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

HAPMAP3_SNPLIST_FOR_GENOMIC_SEM = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="snplist_reference_data",
        sub_group="hapmap",
        sub_folder=PurePath("raw"),
        id=AssetId("hapmap3_snp_list"),
        extension=".snplist",
        filename="w_hm3",
    ),
    url="https://www.dropbox.com/scl/fi/8dsoz94erkkv3fdln1mhb/w_hm3.snplist?rlkey=59lhbaz7bvj3bv669f87fsufg&dl=1",
    md5_hash="e1372a59749eb1f92f7f6931c075f5ac",
)
