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
    url="https://vu.data.surf.nl/public.php/dav/files/zkKbNeNOZAhFXZB/?accept=zip",
    md5_hash="139171d8b859527b6c9231103c4d695f",
)
