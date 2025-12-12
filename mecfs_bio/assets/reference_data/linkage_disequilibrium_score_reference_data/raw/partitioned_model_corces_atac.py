"""
Partitioned LD scores associated from ATAC Seq of immunological cell types.


From the paper:

 "Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types." Nature genetics 50.4 (2018): 621-629.

Data originally derived from:

Corces, M. R. et al. Lineage-specific and single-cell chromatin accessibility
charts human hematopoiesis and leukemia evolution. Nat. Genet. 48,
1193–1203 (2016).

Data is officially located at : https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays

I re-hosted on dropbox.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

PARTITIONED_MODEL_CORCES_ATAC_LD_SCORES_RAW: DownloadFileTask = DownloadFileTask(
    meta=ReferenceFileMeta(
        asset_id="coreces_atac_partitioned_ld_scores",
        group="linkage_disequilibrium_scores",
        sub_group="LDSCORE_LDSC_SEG",
        sub_folder=PurePath("raw"),
        filename="LDSCORE_LDSC_SEG_ldscores_Corces_ATAC_1000Gv3_ldscores",
        extension=".tar",
    ),
    url="https://www.dropbox.com/scl/fi/val9agrddn58ivw5i1ym4/LDSCORE_LDSC_SEG_ldscores_Corces_ATAC_1000Gv3_ldscores.tar?rlkey=4yizbyor0of49jotv8la1ke56&dl=1",
    md5_hash="e39132c2695830a18ef525622c1b037d",
)
