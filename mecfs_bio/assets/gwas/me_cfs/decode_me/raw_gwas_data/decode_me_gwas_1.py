"""
Task to download summary statistics for GWAS-1 from Decode ME via OSF.

GWAS-1 is the primary GWAS in the DECODE ME paper, with the largest number of cases and controls:
See: https://www.medrxiv.org/content/10.1101/2025.08.06.25333109v1.full-text#:~:text=DecodeME%20GWAS%3A%20Overall%20and%20stratified%20analyses
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.osf_retrieve_task import OSFRetrievalTask

DECODE_ME_GWAS_1_TASK = OSFRetrievalTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("decode_me_gwas_1_raw"),
        trait="ME_CFS",
        project="DecodeME",
        sub_dir="raw",
        project_path=PurePath("DecodeME Summary Statistics") / "gwas_1.regenie.gz",
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator=" ")),
    ),
    osf_project_id="rgqs3",
    md5_hash="eabd3c06ffdeb2ec6382bfa67eed7f37",
)
