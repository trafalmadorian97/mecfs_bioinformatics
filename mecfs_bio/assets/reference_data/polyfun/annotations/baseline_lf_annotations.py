"""Reference assets: the allele-bearing baseline-LF 2.2.UKB annotation parquets
and the single sorted annotation matrix built from them.

The annotations come from the ~30GB baselineLF_v2.2.UKB.polyfun.tar.gz bundle,
which is streamed and never stored: only the per-chromosome .annot.parquet
members (~1GB total, carrying A1/A2) are kept. The derived ~19.5M x 187 matrix
and the members directory are path_remap candidates (large, few-file, rarely
read) for relocation to /mnt/d.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.task.annotation_weights.build_baseline_lf_annotation_parquet_task import (
    BuildBaselineLFAnnotationParquetTask,
)
from mecfs_bio.build_system.task.annotation_weights.stream_extract_annotation_parquets_task import (
    StreamExtractAnnotationParquetsTask,
)

# Stream the 30GB polyfun bundle and keep only the per-chromosome .annot.parquet
# members (~1GB, with A1/A2), discarding the ~29GB of LD-score files. The full
# tarball is never stored. This allele-bearing source replaces the allele-less
# .annot.gz path, which forced a lossy per-allele SNP dedup.
BASELINE_LF_ANNOTATION_PARQUET_MEMBERS = StreamExtractAnnotationParquetsTask(
    meta=ReferenceDataDirectoryMeta(
        group="polyfun",
        sub_group="annotations",
        sub_folder=PurePath("raw"),
        id=AssetId("baseline_lf_2.2_ukb_annot_parquet_members"),
    ),
    url="https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/baselineLF_v2.2.UKB.polyfun.tar.gz",
)

BASELINE_LF_ANNOTATION_MATRIX = BuildBaselineLFAnnotationParquetTask.create(
    asset_id="baseline_lf_2.2_ukb_annotations",
    annot_members_task=BASELINE_LF_ANNOTATION_PARQUET_MEMBERS,
)
