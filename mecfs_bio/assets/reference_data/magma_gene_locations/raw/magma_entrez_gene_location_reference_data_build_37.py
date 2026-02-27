from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="magma_reference_data",
        sub_group="gene_locations",
        sub_folder=PurePath("raw"),
        id="magma_gene_location_reference_data_build_37_raw",
        extension=".zip",
        filename="NCBI37.3",
    ),
    url="https://vu.data.surf.nl/public.php/dav/files/Pj2orwuF2JYyKxq/?accept=zip",
    md5_hash=None,
)
