"""Reference assets: the baseline-LF 2.2.UKB annotation tarball and the single
sorted annotation parquet built from it.

The tarball is ~11GB; it and the derived ~19M x 191 parquet are path_remap
candidates (large, few-file, rarely read) for relocation to /mnt/d.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.annotation_weights.build_baseline_lf_annotation_parquet_task import (
    BuildBaselineLFAnnotationParquetTask,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

BASELINE_LF_ANNOTATION_TARBALL = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="polyfun",
        sub_group="annotations",
        sub_folder=PurePath("raw"),
        id=AssetId("baseline_lf_2.2_ukb_annotations_tarball"),
        extension=".tar.gz",
    ),
    url="https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/baselineLF_v2.2.UKB.tar.gz",
    md5_hash=None,  # pin after the first real download (see Step 3)
)

BASELINE_LF_ANNOTATION_MATRIX = BuildBaselineLFAnnotationParquetTask.create(
    asset_id="baseline_lf_2.2_ukb_annotations",
    tarball_task=BASELINE_LF_ANNOTATION_TARBALL,
)
