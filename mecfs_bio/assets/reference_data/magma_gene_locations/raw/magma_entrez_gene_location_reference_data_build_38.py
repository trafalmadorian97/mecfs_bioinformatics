from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

MAGMA_GENE_LOCATION_REFERENCE_DATA_BUILD_38_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="magma_reference_data",
        sub_group="gene_locations",
        sub_folder=PurePath("raw"),
        id="magma_gene_location_reference_data_build_38_raw",
        extension=".zip",
        filename="NCBI38",
    ),
    url="https://vu.data.surf.nl/public.php/dav/files/yj952iHqy5anYhH/?accept=zip",
    md5_hash=None,
)
