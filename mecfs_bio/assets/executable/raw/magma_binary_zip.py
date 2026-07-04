"""
Task to download the MAGMA binary.

The main citation for MAGMA is
 "MAGMA: generalized gene-set analysis of GWAS data."
PLoS computational biology 11.4 (2015): e1004219.

Previously, this Task downloaded MAGMA from the official site.  I changed it to download from my dropbox, since the official
link seemed to have become flaky.

"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.executable_meta import ExecutableMeta
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

MAGMA_1_1_BINARY_ZIPPED = DownloadFileTask(
    meta=ExecutableMeta.create(
        group="gene_set_analysis",
        sub_folder=PurePath("raw"),
        asset_id="magma_1_1_binary_zip",
        extension=".zip",
        filename="magma_v1.10",
    ),
    url="https://www.dropbox.com/scl/fi/yiqx2e3rz1xp79p4882hw/magma_v1.10.zip?rlkey=7nok9rdi6igm74g6t41i8tflv&dl=1",
    md5_hash="139171d8b859527b6c9231103c4d695f",
)
