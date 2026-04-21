from pathlib import PurePath

from mecfs_bio.assets.reference_data.gene_set_data.for_magma.from_gsea_msigdb.full_msigdb_sqlite_zipped import \
    MSGIDB_SQLLITE_RAW
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import ReferenceFileMeta
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
from mecfs_bio.build_system.task.extraction_one_file_from_zip_task import ExtractFromZipTask

MSIGDB_SQLLITE_EXTRACTED = ExtractFromZipTask.create_from_zipped_reference_file(
    source_task=MSGIDB_SQLLITE_RAW,
    asset_id="msigdb_full_sqlite_extracted",
    file_to_extract="msigdb_v2026.1.Hs.db",
    extension=""
)


