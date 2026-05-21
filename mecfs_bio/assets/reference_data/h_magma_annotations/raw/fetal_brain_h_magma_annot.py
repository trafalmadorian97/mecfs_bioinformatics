"""
H-MAGMA gene annotation file for fetal brain tissue.

Maps SNPs to genes using Hi-C chromatin interaction data from fetal brain.
Source: https://github.com/thewonlab/H-MAGMA/tree/master/Input_Files
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

FETAL_BRAIN_H_MAGMA_ANNOT_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="h_magma_annotations",
        sub_group="fetal_brain",
        sub_folder=PurePath("raw"),
        id="fetal_brain_h_magma_annot_raw",
        extension=".annot",
        filename="Fetal_brain.genes",
    ),
    url="https://raw.githubusercontent.com/thewonlab/H-MAGMA/master/Input_Files/Fetal_brain.genes.annot",
    md5_hash="c09bec02e453865eb4444a93c3b030c2",
)
