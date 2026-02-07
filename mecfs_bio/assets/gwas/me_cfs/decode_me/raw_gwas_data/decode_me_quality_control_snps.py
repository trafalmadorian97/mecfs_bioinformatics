"""
Task to get the list of genetic variants passing quality control in DECODE ME.

See the DECODE ME Summary Stats README for an explanation of quality control steps:
https://osf.io/rgqs3/files/axp4k
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.osf_retrieve_task import OSFRetrievalTask

DECODE_ME_QC_SNPS = OSFRetrievalTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("decode_me_snps_passing_quality_control"),
        trait="ME_CFS",
        project="DecodeME",
        sub_dir="raw",
        project_path=PurePath("DecodeME Summary Statistics") / "gwas_qced.var.gz",
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator=" ")),
    ),
    osf_project_id="rgqs3",
)
