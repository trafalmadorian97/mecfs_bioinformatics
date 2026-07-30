"""
H-MAGMA gene annotation file for iPSC-derived astrocytes.

Maps SNPs to genes using Hi-C chromatin interaction data from iPSC-derived
astrocytes.
Source: https://github.com/thewonlab/H-MAGMA/tree/master/Input_Files
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

IPSC_DERIVED_ASTRO_H_MAGMA_ANNOT_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="h_magma_annotations",
        sub_group="ipsc_derived_astro",
        sub_folder=PurePath("raw"),
        id="ipsc_derived_astro_h_magma_annot_raw",
        extension=".annot",
        filename="iPSC_derived_astro.genes",
    ),
    url="https://raw.githubusercontent.com/thewonlab/H-MAGMA/master/Input_Files/iPSC_derived_astro.genes.annot",
    md5_hash="dcef6dd439ac27ebb350c21413266014",
)
