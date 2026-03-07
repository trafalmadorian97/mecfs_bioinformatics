from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import ReferenceFileMeta
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

LAVA_G1000_EUR_REF_LD=DownloadFileTask(
    meta=ReferenceFileMeta(
        group="lava_reference_ld",
        sub_group="g1000_eur",
        sub_folder=PurePath("raw"),
        id="lava_1000g_eur_reference_ld",
        filename="g1000_eur",
        extension=".zip",

    ),
    url="https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download",
    md5_hash=None,
)