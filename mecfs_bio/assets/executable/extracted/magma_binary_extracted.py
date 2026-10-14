from mecfs_bio.assets.executable.raw.magma_binary_zip import MAGMA_1_1_BINARY_ZIPPED
from mecfs_bio.build_system.task.extraction_one_file_from_zip_task import (
    ExtractFromZipTask,
)

MAGMA_1_1_BINARY_EXTRACTED = ExtractFromZipTask.create_from_zipped_executable(
    source_task=MAGMA_1_1_BINARY_ZIPPED,
    asset_id="magma_1_1_binary_extracted",
    file_to_extract="magma",
)
