from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.raw.raw_bellenguez_data import (
    BELLENGUEZ_ET_AL_ALZHIEMERS_RAW,
)
from mecfs_bio.build_system.task.extract_gzip_task import ExtractGzipTextFileTask

BELLINGUEZ_ET_AL_ALZHIEMERS_EXTRACTED = ExtractGzipTextFileTask.create_for_gwas_file(
    source_file_task=BELLENGUEZ_ET_AL_ALZHIEMERS_RAW,
    asset_id="bellinguez_alz_extracted",
)
