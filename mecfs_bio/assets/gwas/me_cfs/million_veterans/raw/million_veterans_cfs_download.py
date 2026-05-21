"""
Task to download the Million Veterans Program CFS summary statistics from GWAS catalog.


See: https://www.ebi.ac.uk/gwas/studies/GCST90479178

"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

MILLION_VETERANS_CFS_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("million_veterans_cfs_eur_raw"),
        trait="ME_CFS",
        project="million_veterans",
        sub_dir="raw",
        project_path=PurePath("GCST90479178.h.tsv.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479178/harmonised/GCST90479178.h.tsv.gz",
    md5_hash="1702d9a7c5caf9488e44f67c776b1ad1",
)
