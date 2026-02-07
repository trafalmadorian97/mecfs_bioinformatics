from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

KEATON_ET_AL_DBP_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("keaton_et_al_dbp_raw"),
        trait="diastolic_blood_pressure",
        project="keaton_et_al_dbp",
        sub_dir="raw",
        project_path=PurePath("GCST90310295.h.tsv.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90310001-GCST90311000/GCST90310295/harmonised/GCST90310295.h.tsv.gz",
    md5_hash="392e09b40dcd4716efb14446745dbb2c",
)
