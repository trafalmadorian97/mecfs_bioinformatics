from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

AEGISDOTTIR_ET_AL_RAW_SYNCOPE_DATA = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("aegisdottir_et_al_raw_syncope_data"),
        trait="syncope",
        project="aegisdottir",
        sub_dir="",
        project_path=PurePath("syncope2023.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://www.dropbox.com/scl/fi/sy1pitlnd1ywe9rsr8qqa/syncope2023.gz?rlkey=crno0suvqolsrp8fsrsjvc17j&dl=1",
    md5_hash="ffe32aefc6f91ad2be0c23a4a1ecbac6",
)
