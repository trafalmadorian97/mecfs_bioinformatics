"""
This is a dataset of partitioned LD scores associated with brain cell types.

It comes from the  paper

 "Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types." Nature genetics 50.4 (2018): 621-629.

 The data is originally derived from transcriptome data from experiments on mouse brains reported in

Cahoy, John D., et al. "A transcriptome database for astrocytes, neurons, and oligodendrocytes: a new resource for understanding brain development and function." Journal of Neuroscience 28.1 (2008): 264-278.

Data is officially located at : https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays

I re-hosted on dropbox.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

PARTITIONED_MODEL_CAHOY_LD_SCORES_RAW: DownloadFileTask = DownloadFileTask(
    meta=ReferenceFileMeta(
        asset_id="cahoy_partitioned_ld_scores",
        group="linkage_disequilibrium_scores",
        sub_group="LDSCORE_LDSC_SEG",
        sub_folder=PurePath("raw"),
        filename="LDSCORE_LDSC_SEG_ldscores_Cahoy_1000Gv3_ldscores",
        extension=".tar",
    ),
    url="https://www.dropbox.com/scl/fi/jnxzjt6wytng8xtdwkgeb/LDSCORE_LDSC_SEG_ldscores_Cahoy_1000Gv3_ldscores.tar?rlkey=ed12zeugsc6q8ddawcy1cifqm&dl=1",
    md5_hash="15c5e5e8f056a12f1da65cd5d45b1963",
)
