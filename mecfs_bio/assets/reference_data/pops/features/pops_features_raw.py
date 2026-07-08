"""
Task to download the full POPs gene-feature collection.

This is the raw combined feature file (PoPS.features.txt.gz) used to generate the
manuscript results. It is NOT in POPs' munged format: it must still be processed by
PopsMungeFeatureDirectoryTask before pops.py can consume it. Hosted by the authors
on Dropbox (linked from the POPs README).

md5_hash is left as None because the file is multiple gigabytes and the authors do
not publish a checksum; pin one here if the download proves unstable.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

POPS_FEATURES_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="pops",
        sub_group="features",
        sub_folder=PurePath("raw"),
        extension=".txt.gz",
        filename="PoPS.features",
        id="pops_features_raw",
    ),
    url=(
        "https://www.dropbox.com/scl/fo/ne7xhxkt4dwhvd52a59ub/"
        "AFKkJu7ACaun1uuE99kmTkc/data/PoPS.features.txt.gz"
        "?rlkey=ltdbcld1enyr1zefg1lfqm61i&dl=1"
    ),
    md5_hash=None,
)
