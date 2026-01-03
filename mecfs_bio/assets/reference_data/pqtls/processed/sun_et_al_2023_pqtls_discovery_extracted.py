"""
Task to extract supplementary table 9 from the Sun et al paper.  This contains the discovery pQTLS
"""

from mecfs_bio.assets.reference_data.pqtls.raw.sun_et_al_2023_pqtls import (
    SUN_ET_AL_2023_PQTL_SUPPLEMENTARY_RAW,
)
from mecfs_bio.build_system.task.extract_sheet_from_excel_file_task import (
    ExtractSheetFromExelFileTask,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat

SUN_ET_AL_2023_DICOVERY_PQTLS_EXTRACTED = ExtractSheetFromExelFileTask.create(
    asset_id="sun_et_al_2023_discovery_pqtls_extracted",
    source_excel_file_task=SUN_ET_AL_2023_PQTL_SUPPLEMENTARY_RAW,
    sheet_name="ST9",  # Table ST9 contains discovery pQTLs,
    out_format=ParquetOutFormat(),
    skiprows=4,
    col_type_mapping={"FitCons_score": str},
)
