from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.raw.raw_johnston_data import (
    JOHNSTON_ET_AL_PAIN_RAW,
)
from mecfs_bio.build_system.task.extract_gzip_task import ExtractGzipTextFileTask

JOHNSTON_ET_AL_PAIN_EXTRACTED = ExtractGzipTextFileTask.create_for_gwas_file(
    source_file_task=JOHNSTON_ET_AL_PAIN_RAW,
    asset_id="johnston_et_al_pain_extracted",

)
