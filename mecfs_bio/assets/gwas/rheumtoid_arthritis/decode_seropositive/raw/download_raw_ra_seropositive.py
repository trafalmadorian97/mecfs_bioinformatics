from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameWhiteSpaceSepTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

DECODE_SEROPOSITIVE_RA_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("decode_ra_seropositive_raw"),
        trait="rheumatoid_arthritis",
        project="decode_ra_seropositive",
        sub_dir="raw",
        project_path=PurePath("RA_Seropositive.txt.gz"),
        read_spec=DataFrameReadSpec(
            DataFrameWhiteSpaceSepTextFormat(
                comment_code="#",
            )
        ),
    ),
    # url="https://www.dropbox.com/scl/fi/8uumw9pwphkyjl2jru03u/RA_Seronegative.txt.gz?rlkey=ad35vhuaeo6oouqsam2nrcxvn&dl=1",
    url="https://www.dropbox.com/scl/fi/d5vxmi04gfiwx8m52xkgh/RA_Seropositive.txt.gz?rlkey=jbo7fvwf27yfdaoodttqo7rld&dl=1",
    # md5_hash="b5c7c41cbd85be11dbc63f5f0402e40f",
    md5_hash=None,
)
