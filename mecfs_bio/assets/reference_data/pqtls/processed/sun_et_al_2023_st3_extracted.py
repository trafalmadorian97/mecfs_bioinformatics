"""
Task to extract supplementary table 3 from the Sun et al paper.

ST3 describes every protein on the Olink Explore 3k platform, including the hg38 coordinates of
its gene, which is what the LDSC tasks use to build cis windows.

The sheet opens with a title row and a blank row before the header, hence skiprows. The
chromosome and coordinate columns are read as strings because the multi-protein complexes on the
platform (see gene_coords.parse_chrom) carry several semicolon-joined values in cells that are
otherwise numeric; keeping them as text lets one parser handle both.
"""

from mecfs_bio.assets.reference_data.pqtls.raw.sun_et_al_2023_pqtls import (
    SUN_ET_AL_2023_PQTL_SUPPLEMENTARY_RAW,
)
from mecfs_bio.build_system.task.dataframe_output import (
    ParquetOutFormat,
)
from mecfs_bio.build_system.task.extract_sheet_from_excel_file_task import (
    ExtractSheetFromExelFileTask,
)
from mecfs_bio.build_system.task.ppp_ldsc.gene_coords import (
    ST3_GENE_CHROM_COL,
    ST3_GENE_END_COL,
    ST3_GENE_START_COL,
)

SUN_ET_AL_2023_ST3_EXTRACTED = ExtractSheetFromExelFileTask.create(
    asset_id="sun_et_al_2023_st3_extracted",
    source_excel_file_task=SUN_ET_AL_2023_PQTL_SUPPLEMENTARY_RAW,
    sheet_name="ST3",  # Table ST3 describes the Olink Explore 3k proteins
    out_format=ParquetOutFormat(),
    skiprows=2,
    col_type_mapping={
        ST3_GENE_CHROM_COL: str,
        ST3_GENE_START_COL: str,
        ST3_GENE_END_COL: str,
    },
)
