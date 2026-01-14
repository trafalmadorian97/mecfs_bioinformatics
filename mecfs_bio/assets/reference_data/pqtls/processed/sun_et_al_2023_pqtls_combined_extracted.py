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
from mecfs_bio.build_system.task.pipes.str_split_exact_col import SplitExactColPipe

SUN_ET_AL_2023_COMBINED_PQTLS_EXTRACTED = ExtractSheetFromExelFileTask.create(
    asset_id="sun_et_al_2023_combined_pqtls_extracted",
    source_excel_file_task=SUN_ET_AL_2023_PQTL_SUPPLEMENTARY_RAW,
    sheet_name="ST10",  # Table ST10 contains combined pQTLs,
    out_format=ParquetOutFormat(),
    skiprows=4,
    col_type_mapping={"FitCons_score": str},
    post_pipe=SplitExactColPipe(
        col_to_split="Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)",
        split_by=":",
        new_col_names=(
            "CHROM_hg37",
            "GENPOS_hg37",
            "A0",
            "A1",
            "imp",
            "v1"
        )
    )
)
