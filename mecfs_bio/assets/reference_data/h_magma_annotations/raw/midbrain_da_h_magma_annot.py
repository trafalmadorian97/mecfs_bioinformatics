"""
H-MAGMA gene annotation file for midbrain dopaminergic neurons.

Maps SNPs to genes using Hi-C chromatin interaction data from midbrain
dopaminergic neurons.
Source: https://github.com/thewonlab/H-MAGMA/tree/master/Input_Files
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

MIDBRAIN_DA_H_MAGMA_ANNOT_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="h_magma_annotations",
        sub_group="midbrain_da",
        sub_folder=PurePath("raw"),
        id="midbrain_da_h_magma_annot_raw",
        extension=".annot",
        filename="Midbrain_DA.genes",
    ),
    url="https://raw.githubusercontent.com/thewonlab/H-MAGMA/master/Input_Files/Midbrain_DA.genes.annot",
    md5_hash="8331ea0cb49248df3f58504400e09aa3",
)
