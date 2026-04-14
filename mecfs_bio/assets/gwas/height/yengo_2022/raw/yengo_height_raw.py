"""
Task to pull raw data from Yengo et al.'s 2022 GWAS of adult height.
Citation:
Yengo, Loïc, et al. "A saturated map of common genetic variants associated
 with human height." Nature 610.7933 (2022): 704-712.

Giant consortium data page: https://giant-consortium.web.broadinstitute.org/index.php/GIANT_consortium_data_files



"""

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
    md5_hash="710ef5ac573e9023b56ec37b959a4e29",
)
