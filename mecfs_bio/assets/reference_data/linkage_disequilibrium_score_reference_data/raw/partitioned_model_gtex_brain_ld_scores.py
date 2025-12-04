"""
Paritioned LD scores associated from brain tissue types.


Comes from the paper:

 "Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types." Nature genetics 50.4 (2018): 621-629.

 Originally derived from the GTEx project

Data is officially located at : https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays

I re-hosted on dropbox.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

PARTITIONED_MODEL_GTEX_BRAIN_LD_SCORES_RAW: DownloadFileTask = DownloadFileTask(
    meta=ReferenceFileMeta(
        asset_id="gtex_brain_partitioned_ld_scores",
        group="linkage_disequilibrium_scores",
        sub_group="LDSCORE_LDSC_SEG",
        sub_folder=PurePath("raw"),
        filename="LDSCORE_LDSC_SEG_ldscores_GTEx_brain_1000Gv3_ldscores",
        extension=".tar",
    ),
    url="https://www.dropbox.com/scl/fi/yg0bf3tq077ut74040wed/LDSCORE_LDSC_SEG_ldscores_GTEx_brain_1000Gv3_ldscores.tar?rlkey=dz8me7skwswwfsybytbe4wg8o&dl=1",
    md5_hash="b9842d0ef3b0b9afe0434d63c43067b2",
)
