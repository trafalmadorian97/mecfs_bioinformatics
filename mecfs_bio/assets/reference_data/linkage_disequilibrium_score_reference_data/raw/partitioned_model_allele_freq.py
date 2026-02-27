"""
This is allele frequency reference data for use with S-LDSC.
See: https://github.com/bulik/ldsc/wiki/Partitioned-Heritability
and
Finucane, Hilary K., et al. "Partitioning heritability by functional annotation using genome-wide association summary statistics." Nature genetics 47.11 (2015): 1228-1235.


Data is officially located at : https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays

"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

LD_SCORE_REGRESSION_ALLELE_FREQS = DownloadFileTask(
    meta=ReferenceFileMeta(
        id="s_ldsc_allele_freqs",
        group="linkage_disequilibrium_score_regression_aux_data",
        sub_group="LDSCORE_1000G_Phase",
        sub_folder=PurePath("raw"),
        filename="LDSCORE_1000G_Phase3_frq",
        extension="tar",
    ),
    url="https://www.dropbox.com/scl/fi/wx987fuu6s3utt6po993e/LDSCORE_1000G_Phase3_frq.tar?rlkey=li2v633jnktnfnkoy00awkzf2&dl=1",
    md5_hash="ac29686ffd5b6378789857a522ebca77",
)
