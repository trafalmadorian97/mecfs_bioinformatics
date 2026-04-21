from mecfs_bio.assets.reference_data.gene_set_data.for_magma.from_gsea_msigdb.raw.full_msigdb_sqlite_zipped import (
    MSGIDB_SQLLITE_RAW,
)
from mecfs_bio.build_system.task.extraction_one_file_from_zip_task import (
    ExtractFromZipTask,
)

MSIGDB_SQLLITE_EXTRACTED = ExtractFromZipTask.create_from_zipped_reference_file(
    source_task=MSGIDB_SQLLITE_RAW,
    asset_id="msigdb_full_sqlite_extracted",
    file_to_extract="msigdb_v2026.1.Hs.db",
    extension="",
)
