"""
H-MAGMA gene annotation file for adult brain tissue.

Maps SNPs to genes using Hi-C chromatin interaction data from adult brain.
Source: https://github.com/thewonlab/H-MAGMA/tree/master/Input_Files
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

ADULT_BRAIN_H_MAGMA_ANNOT_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="h_magma_annotations",
        sub_group="adult_brain",
        sub_folder=PurePath("raw"),
        id="adult_brain_h_magma_annot_raw",
        extension=".annot",
        filename="Adult_brain.genes",
    ),
    url="https://raw.githubusercontent.com/thewonlab/H-MAGMA/master/Input_Files/Adult_brain.genes.annot",
    md5_hash="ee95be1e8e2d5ddb059584559b3bda93",
)
