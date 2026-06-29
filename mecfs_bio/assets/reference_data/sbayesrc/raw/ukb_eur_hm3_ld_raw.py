"""
SBayesRC UKB EUR HapMap3 eigen-decomposed LD reference (raw download).

This is the HapMap3 EUR LD reference from the SBayesRC resource page, used by both
SBayesRC and polypwas.  The archive extracts to a ukbEUR_HM3 directory containing
snp.info, ldm.info and block*.eigen.bin.

Source: https://github.com/zhilizheng/SBayesRC#resources
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

SBAYESRC_UKB_EUR_HM3_LD_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        id="sbayesrc_ukb_eur_hm3_ld_zip",
        group="sbayesrc",
        sub_group="ukb_eur_hm3_ld",
        sub_folder=PurePath("raw"),
        extension=".zip",
    ),
    url="https://gctbhub.cloud.edu.au/data/SBayesRC/resources/v2.0/LD/HapMap3/ukbEUR_HM3.zip",
    md5_hash="773d1001afe8a113b90016a27858b0a2",
)
