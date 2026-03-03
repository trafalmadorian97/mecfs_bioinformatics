from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec, DataFrameTextFormat
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

MILLION_VETERAN_ANKYLOSING_SPONDYLITIS_EUR_DATA_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("million_veterans_spond_eur_raw"),
        trait="ankylosing_spondylitis",
        project="million_veterans",
        sub_dir="raw",
        project_path=PurePath("GCST90476232.h.tsv.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476232/harmonised/GCST90476232.h.tsv.gz",
    md5_hash="4f9c91b7e409aa7b4d92ed9a2de47ced",
)
