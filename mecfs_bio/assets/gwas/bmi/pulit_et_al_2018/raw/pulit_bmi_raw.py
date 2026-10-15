"""
Task to download the GIANT + UK Biobank meta-analysis GWAS of body mass index (sexes combined) from:
"Meta-analysis of genome-wide association studies for body fat distribution in 694,649 individuals of European ancestry." Pulit, SL et al. Human Molecular Genetics, 2019.
Files are hosted at: https://zenodo.org/records/1251813
"""

from pathlib import PurePath

import polars as pl

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

PULIT_ET_AL_2018_BMI_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("pulit_et_al_2018_bmi_raw"),
        trait="bmi",
        project="pulit_et_al_2018",
        sub_dir="raw",
        read_spec=DataFrameReadSpec(
            format=DataFrameTextFormat(
                separator=" ", null_values=["NA"], schema_overrides={"N": pl.Float64()}
            )
        ),
        project_path=PurePath("bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz"),
    ),
    url="https://zenodo.org/records/1251813/files/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1",
    md5_hash="9c0d6c5249ffa816681a92f8c4bba424",
)
