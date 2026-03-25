from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

XUE_WHOLE_BRAINSTEM_VOLUME_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("xue_et_al_whole_brainstem_raw"),
        trait="whole_brainstem_volume",
        project="xue_et_al_2025",
        sub_dir="raw",
        project_path=PurePath("eur_whole_brainstem_volume.gz"),
        read_spec=DataFrameReadSpec(DataFrameTextFormat(separator="\t")),
    ),
    url="https://zenodo.org/records/13382122/files/eur_whole_brainstem_volume.gz?download=1",
    md5_hash="db15cd75b74e03de93b13711d6ad1e85",
)
