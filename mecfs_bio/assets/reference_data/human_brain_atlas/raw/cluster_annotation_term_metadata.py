"""
Task to fetch some metadata from the human brain atlas.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

CLUSTER_ANNOTATION_TERM_METADATA = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="human_brain_atlas",
        sub_group="metadata",
        sub_folder=PurePath("raw"),
        id=AssetId("human_brain_atlas_cluster_annotation_Term_metadata"),
        filename="cluster_annotation_term",
        extension=".csv",
        read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
    ),
    url="https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/metadata/WHB-taxonomy/20240330/cluster_annotation_term.csv",
    md5_hash="2533ac1498b35a1ce386a3d97137e9c3",
)
