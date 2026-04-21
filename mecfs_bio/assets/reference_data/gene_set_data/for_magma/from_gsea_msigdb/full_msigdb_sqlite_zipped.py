from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import ReferenceFileMeta
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

MSGIDB_SQLLITE_RAW=DownloadFileTask(
    meta=ReferenceFileMeta(
    id = AssetId("msigdb_full_sqlite"),
    group = "gene_set_data",
    sub_group = "msigdb",
    sub_folder = PurePath("raw"),
    filename = "msigdb.v2026.1.Hs",
    extension = ".zip",
),
    url="https://www.dropbox.com/scl/fi/oosi9pe56ip7f9q3jau8f/msigdb_v2026.1.Hs.db.zip?rlkey=xwd366dg1suwgab8bkzs6hlf2&dl=1",
    md5_hash=None

)