from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

PARTITIONED_MODEL_SUPERCLUSTERS_YAO: DownloadFileTask = DownloadFileTask(
    meta=ReferenceFileMeta(
        asset_id="yao_superclusters_partitioned_ld_scores",
        group="linkage_disequilibrium_scores",
        sub_group="yao_partitioned_ldscores",
        sub_folder=PurePath("raw"),
        filename="superclusters",
        extension=".tar",
    ),
    url="https://zenodo.org/records/10628706/files/superclusters.tar?download=1",
    md5_hash=None,
)
