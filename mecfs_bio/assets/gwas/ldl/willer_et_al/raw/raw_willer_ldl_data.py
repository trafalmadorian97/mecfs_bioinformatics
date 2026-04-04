from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec, DataFrameTextFormat
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

WILLER_LDL_EUR_DATA_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("willer_ldl_eur_raw"),
        trait="ldl",
        project="willer_et_al_2013",
        sub_dir="raw",
        project_path=PurePath("24097068-GCST002222-EFO_0004611.h.tsv.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t",
                                                               null_values=["NA"]
                                                               )),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST002001-GCST003000/GCST002222/harmonised/24097068-GCST002222-EFO_0004611.h.tsv.gz",
    md5_hash="260072772c7a1c9ddebad23ff3a74677",

)
