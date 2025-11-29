"""
This is one of the reference datasets from the Stratified Linkage Disequilibrium Score Regression Paper.
It contains SNP-weights for use in the regression.

These weights attempt to adjust for two issues:
- The chi squared statistics used by S-LDSC as regression targets are not independent
- Heteroskedasticity: chi-squared statistics corresponding to SNPs with higher LD scores will have higher variance

see:
Finucane, Hilary K., et al. "Partitioning heritability by functional annotation using genome-wide association summary statistics." Nature genetics 47.11 (2015): 1228-1235.

https://github.com/bulik/ldsc/wiki/Partitioned-Heritability




Data is officially located at : https://console.cloud.google.com/storage/browser/broad-alkesgroup-public-requester-pays
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

HM3_LD_SCORE_REGRESSION_WEIGHTS = DownloadFileTask(
    meta=ReferenceFileMeta(
        asset_id="s_ldsc_hapmap3_regression_weights",
        group="linkage_disequilibrium_score_regression_aux_data",
        sub_group="LDSCORE_1000G_Phase",
        sub_folder=PurePath("raw"),
        filename="LDSCORE_1000G_Phase3_weights_hm3_no_MHC",
        extension=".tar",
    ),
    url="https://www.dropbox.com/scl/fi/wycm3egtm4if69wcs65pb/LDSCORE_1000G_Phase3_weights_hm3_no_MHC.tar?rlkey=e11bcg1rpkh4d7subfwjixfpj&dl=1",
    md5_hash="a98ac0f089ee285177544a3e6e721ca3",
)
