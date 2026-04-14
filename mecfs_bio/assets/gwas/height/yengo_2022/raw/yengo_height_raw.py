from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

YENGO_ET_AL_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("yengo_height_raw"),
        trait="height",
        project="yengo_et_al_2022",
        sub_dir="raw",
        project_path=PurePath("GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz"),
        read_spec=DataFrameReadSpec(DataFrameTextFormat(separator="\t")),
    ),
    url="https://giant-consortium.web.broadinstitute.org/images/f/f7/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz",
    md5_hash=None,
)
