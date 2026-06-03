"""
Version 1.1 of the 'baseline' partitioned-model LD scores from the S-LDSC authors.

This is an older baseline model than the v1.2 file used by default (see
`partitioned_model_baseline_ld_scores.py`).  It is included to investigate how much the
baseline model version contributes to differences in S-LDSC results between our pipeline and
collaborators following the standard LDSC cell-type tutorial, which historically points at an
unversioned `1000G_Phase3_baseline` directory.

Since this dataset is hosted in a requester-pays Google Cloud bucket, it has been mirrored in
Dropbox.

Data is officially located at : https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

THOUSAND_GENOME_BASE_MODEL_PARTITIONED_LD_SCORES_V1_1 = DownloadFileTask(
    meta=ReferenceFileMeta(
        id="thousand_genomes_baseline_model_partitioned_ld_scores_v1_1",
        group="linkage_disequilibrium_scores",
        sub_group="LDSCORE_1000G_Phase3_baseline_v1_1",
        sub_folder=PurePath("raw"),
        filename="LDSCORE_1000G_Phase3_baseline_v1.1_ldscores",
        extension=".tar",
    ),
    url="https://www.dropbox.com/scl/fi/jq22197f51pyxvtogdu5v/LDSCORE_1000G_Phase3_baseline_v1.1_ldscores.tar?rlkey=tn2ef127jz5qjcbs4tp3dtgbq&dl=1",
    md5_hash="02101555c13f760a4f0608905f569cf8",
)
