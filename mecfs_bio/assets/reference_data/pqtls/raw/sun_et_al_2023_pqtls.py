"""
Access the supplementary material from the paper Sun, Benjamin B., et al. "Plasma proteomic associations with genetics and health in the UK Biobank." Nature 622.7982 (2023): 329-338.

Contains tables of pQTLs from the UK Biobank Pharma Proteomics  Project


"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

SUN_ET_AL_2023_PQTL_SUPPLEMENTARY_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="pqtl_data",
        sub_group="sun_et_al_2023",
        sub_folder=PurePath("raw"),
        id="sun_et_al_2023_pqtl",
        filename="sun_et_al_2023_supplementary",
        extension=".xlsx",
    ),
    url="https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06592-6/MediaObjects/41586_2023_6592_MOESM3_ESM.xlsx",
    md5_hash="aeb2acdc5586034e91197658a10a7d70",
)
