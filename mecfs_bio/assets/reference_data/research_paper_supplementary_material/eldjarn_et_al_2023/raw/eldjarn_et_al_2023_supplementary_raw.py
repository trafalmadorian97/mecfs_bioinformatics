"""
Task to download the supplementary tables of Eldjarn, Grimur H., et al. "Large-scale plasma
proteomics comparisons through genetics and disease associations." Nature 622.7982 (2023):
348-358, the deCODE SomaScan and Olink pQTL comparison.

ST36 lists the highly pleiotropic pQTLs, used to flag deCODE coloc results. ST29 lists the
SomaScan assays with deCODE summary statistics; the deCODE manifest script reads it from its
own download (see regenerate_decode_manifest.py), since the manifest decides which protein
Tasks exist.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

ELDJARN_ET_AL_2023_SUPPLEMENTARY_URL = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06563-x/MediaObjects/41586_2023_6563_MOESM3_ESM.xlsx"
ELDJARN_ET_AL_2023_SUPPLEMENTARY_MD5 = "9d2d3a94d4b3541488936cf0ca93e093"

ELDJARN_ET_AL_2023_SUPPLEMENTARY_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="research_paper_supplementary_material",
        sub_group="eldjarn_et_al_2023",
        sub_folder=PurePath("raw"),
        id="eldjarn_et_al_2023_supplementary_raw",
        filename="41586_2023_6563_MOESM3_ESM",
        extension=".xlsx",
    ),
    url=ELDJARN_ET_AL_2023_SUPPLEMENTARY_URL,
    md5_hash=ELDJARN_ET_AL_2023_SUPPLEMENTARY_MD5,
)
