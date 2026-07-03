from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

MAGMA_EUR_BUILD_37_1K_GENOMES_REF = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="magma_reference_data",
        sub_group="magma_ld_ref",
        sub_folder=PurePath("raw"),
        id=AssetId("magma_eur_1k_genomes_build_37_ld_ref_raw"),
        extension=".zip",
        filename="g1000_eur",
    ),
    url="https://www.dropbox.com/scl/fi/9t9jdec7nlhll73chdehi/g1000_eur.zip?rlkey=6rdirky37vksdsfxi9d8gl27i&dl=1",
    md5_hash="1919cb5c79bbe7871aed71ae4abe6217",
)
