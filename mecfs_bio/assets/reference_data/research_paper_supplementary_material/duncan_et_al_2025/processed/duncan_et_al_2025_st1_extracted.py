"""
Extract supplementary table 1 from the Duncan et al supplementary material

This table includes cluster weights for schizophrenia, as well as cluster descriptions
I think we may be able to reuse the cluster descriptions for our own analysis.
"""

from mecfs_bio.assets.reference_data.research_paper_supplementary_material.duncan_et_al_2025.raw.duncan_et_al_2025_supplementary_raw import (
    DUNCAN_ET_AL_2025_SUPPLEMENTARY_RAW,
)
from mecfs_bio.build_system.task.extract_sheet_from_excel_file_task import (
    ExtractSheetFromExelFileTask,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat

DUNCAN_ET_AL_2025_ST1_EXTRACTED = ExtractSheetFromExelFileTask.create(
    asset_id="duncan_et_al_2025_st1_extracted",
    source_excel_file_task=DUNCAN_ET_AL_2025_SUPPLEMENTARY_RAW,
    sheet_name="Supp_Table_1_SCHIZOPHRENIA",  # Table ST1 contains cluster weights for schizophrenia, as well as cluster descriptions,
    out_format=ParquetOutFormat(),
    col_type_mapping={
        "Class auto-annotation": str,
        "Neurotransmitter auto-annotation": str,
        "Subtype auto-annotation": str,
        "Neuropeptide auto-annotation": str,
    },
)
