"""
Tasks to download the full set of MSIG_DB Gene sets and convert them to a table parquet.


"""
from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
from mecfs_bio.build_system.task.json_to_parquet_task import JsonToParquetTask
from mecfs_bio.build_system.task.unpack_map_parquet_task import UnpackMapParquetTask

MSIGDB_JSON_GENE_SETS = DownloadFileTask(
    md5_hash="e190568dc54bfff94b6a32537bb0725e",
    url="https://www.dropbox.com/scl/fi/vb6jzppnvtcvmfe3aym7e/msigdb.v2026.1.Hs.txt?rlkey=1457ingl0ku0ns3gib5uttpih&dl=0",
    meta=ReferenceFileMeta(
        id=AssetId("msigdb_full_json"),
        group="gene_set_data",
        sub_group="msigdb",
        sub_folder=PurePath("raw"),
        filename="msigdb.v2026.1.Hs",
        extension=".txt",
    ),
)

MSIGDB_JSON_GENE_SETS_PARQUET = JsonToParquetTask.create(
    json_task=MSIGDB_JSON_GENE_SETS,
    asset_id="msigdb_full_json_parquet",
)

MSIGDB_GENE_SETS_PARQUET_UNPACKED = UnpackMapParquetTask.create(
    source_task=MSIGDB_JSON_GENE_SETS_PARQUET,
    asset_id="msigdb_gene_sets_parquet_unpacked",
)
