"""
H-MAGMA gene annotation file for cortical neurons.

Maps SNPs to genes using Hi-C chromatin interaction data from cortical neurons.
Source: https://github.com/thewonlab/H-MAGMA/tree/master/Input_Files
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

CORTICAL_NEURON_H_MAGMA_ANNOT_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="h_magma_annotations",
        sub_group="cortical_neuron",
        sub_folder=PurePath("raw"),
        id="cortical_neuron_h_magma_annot_raw",
        extension=".annot",
        filename="Cortical_Neuron.genes",
    ),
    url="https://raw.githubusercontent.com/thewonlab/H-MAGMA/master/Input_Files/Cortical_Neuron.genes.annot",
    md5_hash="ad152709626a837ed07d10c94deff806",
)
