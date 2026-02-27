"""
Partitioned LD scores associated from immunological cell types.


From the paper:

 "Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types." Nature genetics 50.4 (2018): 621-629.

Data originally derived from:

Heng, Tracy SP, et al. "The Immunological Genome Project: networks of gene expression in immune cells." Nature immunology 9.10 (2008): 1091-1094.

Data is officially located at : https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays

I re-hosted on dropbox.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

PARTITIONED_MODEL_IMMGEN_LD_SCORES_RAW: DownloadFileTask = DownloadFileTask(
    meta=ReferenceFileMeta(
        id="immgen_partitioned_ld_scores",
        group="linkage_disequilibrium_scores",
        sub_group="LDSCORE_LDSC_SEG",
        sub_folder=PurePath("raw"),
        filename="LDSCORE_LDSC_SEG_ldscores_ImmGen_1000Gv3_ldscores",
        extension=".tar",
    ),
    url="https://www.dropbox.com/scl/fi/kb30br7qcc68lsn84zak0/LDSCORE_LDSC_SEG_ldscores_ImmGen_1000Gv3_ldscores.tar?rlkey=rndsvy5q7m14sxp3swgu4r0lt&dl=1",
    md5_hash="0cddb5ac5348045e4adabe9dac1c9485",
)
