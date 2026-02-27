"""
GWAS summary statistics from the 2023 meta-analysis of Liu et al.
See: https://www.nature.com/articles/s41588-023-01384-0
Google Drive link to downloads: https://drive.google.com/drive/folders/1MbjKuYp77SQEb1Q06nQG2Asq7gt6kGLO
"""

from pathlib import PurePath

import polars as pl

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_from_google_drive_task import (
    DownloadFromGoogleDriveTask,
)

LIU_ET_AL_2023_IBD_META = DownloadFromGoogleDriveTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("liu_et_all_2023_ibd_meta_raw"),
        trait="inflammatory_bowl_disease",
        project="liu_et_al_2023_ibd_meta",
        sub_dir="raw",
        project_path=PurePath("sumstats/ibd_EAS_EUR_SiKJEF_meta_IBD.TBL.txt.gz"),
        read_spec=DataFrameReadSpec(
            DataFrameTextFormat(separator="\t", schema_overrides={"CHR": pl.String()})
        ),
    ),
    file_id="1KeYZd7WIqdQjdbRfUp084JDetN1S_01e",
)
