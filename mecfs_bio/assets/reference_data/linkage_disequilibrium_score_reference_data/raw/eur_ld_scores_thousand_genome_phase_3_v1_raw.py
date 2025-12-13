from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

"""
Note: reference LD scores derived from the 1000 genome project were formerly stored on the broad institute's website, but have since been move.
See: https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/README_new_data_location.txt
New data location: https://zenodo.org/records/7768714
"""
THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        asset_id="thousand_genomes_phase_3_v1_eur_ld_scores_tar_gzip",
        group="linkage_disequilibrium_scores",
        sub_group="thousand_genomes_phase_3_v1",
        sub_folder=PurePath("raw"),
        extension=".tar.gz",
    ),
    url="https://zenodo.org/records/7768714/files/1000G_Phase3_ldscores.tgz?download=1",
    md5_hash="e4352ccf778f296835d73985350a863b",
)
