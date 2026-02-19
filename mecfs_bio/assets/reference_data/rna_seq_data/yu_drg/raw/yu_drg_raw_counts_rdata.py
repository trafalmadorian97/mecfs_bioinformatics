"""
Task to RDATA object containing RNAseq count data from the paper "Leveraging deep single-soma RNA sequencing to explore the neural basis of human somatosensation" by Yu et al.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

YU_DRG_SRC1_RDATA = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="rna_seq_data",
        sub_group="yu_dorsal_root_ganglia",
        sub_folder=PurePath("raw"),
        id="yu_drg_source_1_rdata",
        filename="Hs_LCM_final",
        extension=".RData",
    ),
    url="https://github.com/taimeimiaole/NN_hDRG-neuron-sequencing/raw/refs/heads/main/Source_code_1/Hs_LCM_final.RData",
    md5_hash="3064b44e1ee61c8d7835114adbf75d31",
)
