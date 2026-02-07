from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

"""
Official page:
https://www.gsea-msigdb.org/gsea/downloads.jsp
File is labeled "Human Gene Set GMT
file set (ZIPped)"
Official page does not provide a direct download link, so hosted on Dropbox.
"""
GSEA_HUMAN_GENE_SET_2025_V1 = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="gene_set_data",
        sub_group="gsea_human_data",
        sub_folder=PurePath("raw"),
        id=AssetId("gsea_human_data_2025_v1"),
        filename="msigdb_v2025.1",
        extension=".zip",
    ),
    url="https://www.dropbox.com/scl/fi/4902crqdai9pmleiid4q5/msigdb_v2025.1.Hs_files_to_download_locally.zip?rlkey=12ib2iy5pzwxs3eauli4h1gip&st=iojgsqok&dl=1",  # official page
    md5_hash="eef0cc56521fecb1f0949733d211b8d2",
)
