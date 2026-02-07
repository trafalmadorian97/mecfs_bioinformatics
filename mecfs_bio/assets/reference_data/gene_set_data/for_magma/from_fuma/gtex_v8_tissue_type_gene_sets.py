from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

"""
Official page:
https://fuma.ctglab.nl/downloadPage
Labeled: GTEx v8 specific expression differentially expressed genes

Official page does not provide a direct download link, so I hosted on Dropbox.
"""
GTEx_V8_TISSUE_EXPRESSION_BASED_GENE_SET_DATA = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="gene_set_data",
        sub_group="fuma_gtex_human_data",
        sub_folder=PurePath("raw"),
        id=AssetId("gtex_v8_tissue_expression_data"),
        filename="gtex_v8_ts_DEG",
        extension=".txt",
    ),
    url="https://www.dropbox.com/scl/fi/qpsg0l5nw1sngy04otzcj/gtex_v8_ts_DEG.txt?rlkey=ku0ittmmkehyuldiqy8wvl22w&dl=1",
    md5_hash="9e0ff40fcb2dec7603ff335efc352f24",
)
