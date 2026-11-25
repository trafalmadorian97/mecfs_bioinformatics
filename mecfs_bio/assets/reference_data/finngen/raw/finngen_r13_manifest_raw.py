"""
Task to download the manifest of FinnGen release 13 summary statistics.

The manifest lists one row per endpoint, with its phenocode, description, category,
case and control counts, and the location of its summary statistics file.

Stored alongside the summary statistics in FinnGen's public bucket.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

FINNGEN_R13_MANIFEST_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="finngen",
        sub_group="r13_manifest",
        sub_folder=PurePath("raw"),
        id="finngen_r13_manifest_raw",
        extension=".tsv",
        filename="finngen_R13_manifest",
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv",
    md5_hash="a70f4ecb57acf7c23e8f5d71be11d936",
)
