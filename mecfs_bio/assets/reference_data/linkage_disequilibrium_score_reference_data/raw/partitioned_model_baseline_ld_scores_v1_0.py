"""
Version 1.0 (original, unversioned) of the 'baseline' partitioned-model LD scores from the S-LDSC authors.

This is the oldest baseline model (older than v1.1 and v1.2). The tarball's inner directory is
`1000G_EUR_Phase3_baseline/`, which matches the path in collaborator Martin's LDSC command
(`/LDSC/1000G_EUR_Phase3_baseline/baseline.`) -- so v1.0 is the prime candidate for the baseline
he actually used. Included to test whether matching his baseline moves our cell-type S-LDSC
p-values toward his.

Since this dataset is hosted in a requester-pays Google Cloud bucket, it has been mirrored in
Dropbox.

Data is officially located at : https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

THOUSAND_GENOME_BASE_MODEL_PARTITIONED_LD_SCORES_V1_0 = DownloadFileTask(
    meta=ReferenceFileMeta(
        id="thousand_genomes_baseline_model_partitioned_ld_scores_v1_0",
        group="linkage_disequilibrium_scores",
        sub_group="LDSCORE_1000G_Phase3_baseline_v1_0",
        sub_folder=PurePath("raw"),
        filename="LDSCORE_1000G_Phase3_baseline_ldscores",
        extension=".tar",
    ),
    url="https://www.dropbox.com/scl/fi/s1rt6pg00u1tzvu2hu61s/LDSCORE_1000G_Phase3_baseline_ldscores.tar?rlkey=jo0puelwzrkrcbewiptt1ps61&dl=1",
    md5_hash="b6999fd3cac974a40184416afc12e397",
)
