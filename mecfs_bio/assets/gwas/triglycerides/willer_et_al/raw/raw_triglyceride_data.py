"""
Task to download summary statistics of Willer et al.'s GWAS of triglycerides.

GWAS catalog study page: https://www.ebi.ac.uk/gwas/publications/24097068


Citation: "Discovery and refinement of loci associated with lipid levels." Nature genetics 45, no. 11 (2013): 1274-1283.

GWAS catalog FTP https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST002001-GCST003000/GCST002216
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

WILLER_TG_EUR_DATA_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("willer_tg_eur_raw"),
        trait="triglycerides",
        project="willer_et_al_2013",
        sub_dir="",
        project_path=PurePath("24097068-GCST002216-EFO_0004530.h.tsv.gz"),
        read_spec=DataFrameReadSpec(
            format=DataFrameTextFormat(separator="\t", null_values=["NA"])
        ),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST002001-GCST003000/GCST002216/harmonised/24097068-GCST002216-EFO_0004530.h.tsv.gz",
    md5_hash="4cf3f4cc06796d289a05a5fe5c7030a1",
)
