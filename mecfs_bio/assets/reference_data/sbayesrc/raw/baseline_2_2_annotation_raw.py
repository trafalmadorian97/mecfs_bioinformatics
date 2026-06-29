"""
SBayesRC baseline-LD model 2.2 functional annotation (raw download).

The archive contains a single TAB delimited file annot_baseline2.2.txt with
columns SNP, Intercept and 96 functional annotations for ~8M SNPs (baseline model
2.2; Marquez-Luna 2021, doi:10.1038/s41467-021-25171-9).  This is the optional
annotation input to SBayesRC.

Source: https://github.com/zhilizheng/SBayesRC#resources
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

SBAYESRC_BASELINE_2_2_ANNOTATION_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        id="sbayesrc_baseline_2_2_annotation_zip",
        group="sbayesrc",
        sub_group="baseline_2_2_annotation",
        sub_folder=PurePath("raw"),
        extension=".zip",
    ),
    url="https://gctbhub.cloud.edu.au/data/SBayesRC/resources/v2.0/Annotation/annot_baseline2.2.zip",
    md5_hash=None,
)
