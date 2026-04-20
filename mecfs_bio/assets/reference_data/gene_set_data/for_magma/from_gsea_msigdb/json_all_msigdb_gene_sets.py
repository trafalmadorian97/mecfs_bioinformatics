from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import ReferenceFileMeta
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

MSIGDB_JSON_GENE_SETS = DownloadFileTask(
    md5_hash=None,
    url="https://www.dropbox.com/scl/fi/vb6jzppnvtcvmfe3aym7e/msigdb.v2026.1.Hs.txt?rlkey=1457ingl0ku0ns3gib5uttpih&dl=0",
    meta=ReferenceFileMeta(
        id=AssetId("msigdb_full_json"),
        group = "gene_set_data",
        sub_group="msigdb",
        sub_folder=PurePath("raw"),
        filename="msigdb.v2026.1.Hs",
        extension=".txt",
)

)