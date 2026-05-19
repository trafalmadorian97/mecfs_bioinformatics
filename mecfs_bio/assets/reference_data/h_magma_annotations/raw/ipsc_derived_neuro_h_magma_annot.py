"""
H-MAGMA gene annotation file for iPSC-derived neurons.

Maps SNPs to genes using Hi-C chromatin interaction data from iPSC-derived
neurons.
Source: https://github.com/thewonlab/H-MAGMA/tree/master/Input_Files

Note: the upstream file in the H-MAGMA repository has a leading space in its
name (" iPSC_derived_neuro.genes.annot"), which is reflected in the download
URL as `%20`. The local asset name has no such leading space.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

IPSC_DERIVED_NEURO_H_MAGMA_ANNOT_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="h_magma_annotations",
        sub_group="ipsc_derived_neuro",
        sub_folder=PurePath("raw"),
        id="ipsc_derived_neuro_h_magma_annot_raw",
        extension=".annot",
        filename="iPSC_derived_neuro.genes",
    ),
    url="https://raw.githubusercontent.com/thewonlab/H-MAGMA/master/Input_Files/%20iPSC_derived_neuro.genes.annot",
    md5_hash="40ca1fd7d124d51bfe7ff589afd65b8c",
)
