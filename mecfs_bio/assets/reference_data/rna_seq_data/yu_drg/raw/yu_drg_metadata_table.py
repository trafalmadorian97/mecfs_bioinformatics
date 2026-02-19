"""
Task to download cluster metadata from the paper "Leveraging deep single-soma RNA sequencing to explore the neural basis of human somatosensation" by Yu et al.
Main github repo for this paper is here: https://github.com/taimeimiaole/NN_hDRG-neuron-sequencing/
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

YU_DRG_METADATA_TABLE = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="rna_seq_data",
        sub_group="yu_dorsal_root_ganglia",
        sub_folder=PurePath("raw"),
        id="yu_drg_source_3_metadata_table",
        filename="human_drg_meta_data_newname",
        extension=".csv",
        read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
    ),
    url="https://raw.githubusercontent.com/taimeimiaole/NN_hDRG-neuron-sequencing/refs/heads/main/Source_code_3/human_drg_meta_data_newname.csv",
    md5_hash="db997c9076aca4d38e0646a254670974",
)
