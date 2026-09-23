"""Reference asset: the baseline-LF 2.2.UKB per-chromosome annotation LD-score parquets.

Streamed from the same ~30GB baselineLF_v2.2.UKB.polyfun.tar.gz bundle as the
annotation matrix, keeping the .l2.ldscore.parquet members (each: CHR, SNP, BP,
A1, A2 + 187 annotation-named LD-score columns). The members directory is a
path_remap candidate (~29GB, few files, rarely read).
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.task.annotation_weights.stream_extract_annotation_parquets_task import (
    LDSCORE_PARQUET_MEMBER_RE,
    StreamExtractAnnotationParquetsTask,
)

BASELINE_LF_ANNOTATION_LDSCORE_MEMBERS = StreamExtractAnnotationParquetsTask(
    meta=ReferenceDataDirectoryMeta(
        group="polyfun",
        sub_group="annotations",
        sub_folder=PurePath("raw"),
        id=AssetId("baseline_lf_2.2_ukb_ldscore_parquet_members"),
    ),
    url="https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/baselineLF_v2.2.UKB.polyfun.tar.gz",
    member_pattern=LDSCORE_PARQUET_MEMBER_RE,
    dest_stem="baselineLF2.2.UKB.{chrom}.l2.ldscore",
)
