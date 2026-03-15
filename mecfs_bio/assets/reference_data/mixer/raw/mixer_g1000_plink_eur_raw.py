"""
Task to download European 1000 genomes reference data for use as an LD reference for MiXeR.

"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

MIXER_RAW_G1000_PLINK_DATA = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="mixer",
        sub_group="g1000_plink_reference",
        sub_folder=PurePath("raw"),
        extension=".zip",
        id="mixer_raw_g1000_plink_reference_data",
    ),
    url="https://www.dropbox.com/scl/fi/5pmw6ynentxq408dv5e6j/1000G_EUR_Phase3_plink.zip?rlkey=njkgknz1ntwf7v3nk5uivaxve&dl=1",
    md5_hash=None,  # "07c7ad862cfc612b1049c7fa25a5d16d"
)
