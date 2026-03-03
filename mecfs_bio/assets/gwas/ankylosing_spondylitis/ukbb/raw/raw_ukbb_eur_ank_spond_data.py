"""
Summary statistics from Ankylosing spondylitis in the UK biobank whole genome sequencing study


See: https://www.ebi.ac.uk/gwas/studies/GCST90474065
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

UK_BIOBANK_ANKYLOSING_SPONDYLITIS_EUR_DATA_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("uk_biobank_spond_eur_raw"),
        trait="ankylosing_spondylitis",
        project="uk_biobank",
        sub_dir="raw",
        project_path=PurePath("GCST90474065.h.tsv.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90474001-GCST90475000/GCST90474065/harmonised/GCST90474065.h.tsv.gz",
    md5_hash="de05204abea442b35338963d6b571372",
)
