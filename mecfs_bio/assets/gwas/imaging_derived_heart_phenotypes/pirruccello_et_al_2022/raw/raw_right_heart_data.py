from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

PIRRUCCELLO_RAW_RIGHT_HEART_DATA = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        short_id=AssetId("pirruccello_et_al_2022_raw_right_heart_data"),
        trait="imaging_derived_heart_phenotypes",
        project="pirruccello_et_al",
        sub_dir="raw",
        project_path=PurePath("Pirruccello_2022_UKBB_HeartStructures.zip"),
        # read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://personal.broadinstitute.org/ryank/Pirruccello_2022_UKBB_HeartStructures.zip",
    md5_hash=None,
)
