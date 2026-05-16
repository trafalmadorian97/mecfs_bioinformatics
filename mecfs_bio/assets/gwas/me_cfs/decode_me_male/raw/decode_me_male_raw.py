from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec, DataFrameTextFormat
from mecfs_bio.build_system.task.osf_retrieve_task import OSFRetrievalTask

DECODE_ME_GWAS_MALE_TASK = OSFRetrievalTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("decode_me_gwas_male_raw"),
        trait="ME_CFS",
        project="DecodeME",
        sub_dir="raw",
        project_path=PurePath("DecodeME Summary Statistics") / "gwas_1_male.regenie.gz",
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator=" ")),
    ),
    osf_project_id="rgqs3",
    md5_hash=None,
)
