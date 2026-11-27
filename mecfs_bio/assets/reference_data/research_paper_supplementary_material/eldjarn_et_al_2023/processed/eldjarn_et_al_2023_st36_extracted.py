"""
Extract sheet ST36_pleiotropic_pQTLs (highly pleiotropic pQTLs) from the Eldjarn et al. 2023
supplementary tables, as the sheet reads: title row, blank row, then a header whose names wrap
onto several lines.
"""

from mecfs_bio.assets.reference_data.research_paper_supplementary_material.eldjarn_et_al_2023.raw.eldjarn_et_al_2023_supplementary_raw import (
    ELDJARN_ET_AL_2023_SUPPLEMENTARY_RAW,
)
from mecfs_bio.build_system.task.dataframe_output import ParquetOutFormat
from mecfs_bio.build_system.task.extract_sheet_from_excel_file_task import (
    ExtractSheetFromExelFileTask,
)

ELDJARN_ET_AL_2023_ST36_EXTRACTED = ExtractSheetFromExelFileTask.create(
    asset_id="eldjarn_et_al_2023_st36_extracted",
    source_excel_file_task=ELDJARN_ET_AL_2023_SUPPLEMENTARY_RAW,
    sheet_name="ST36_pleiotropic_pQTLs",
    out_format=ParquetOutFormat(),
    skiprows=2,
)
