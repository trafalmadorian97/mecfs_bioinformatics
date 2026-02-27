"""
Download genetic variant labels for the Broad institute's UK Biobank LD matrix
at Chromosome 1 173000001-17600000 (GRCh37)
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

CHR1_173000001_17600000_UKBB_LD_LABELS_DOWNLOAD = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="ukbb_reference_ld",
        sub_group="chr1_173000001_17600000",
        sub_folder=PurePath("raw"),
        id=AssetId("ukbb_chr1_173000001_17600000_ld_labels"),
        filename="chr1_173000001_176000001",
        extension=".gz",
        read_spec=DataFrameReadSpec(DataFrameTextFormat(separator="\t")),
    ),
    url="https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/chr1_173000001_176000001.gz",
    md5_hash=None,
)
