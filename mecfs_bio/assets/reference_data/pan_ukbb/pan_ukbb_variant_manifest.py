"""
The Pan-UK Biobank variant manifest (full_variant_qc_metrics.txt.bgz, 2.7 GB, GRCh37),
pinned by md5 (equal to its S3 ETag; last modified 2020-08-28).

It stays in the asset store as a dependency of the Pan-UKBB panel; its own sub_folder
(raw/) lets a path_remap rule move it to another disk.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
from mecfs_bio.build_system.task.genome_reference_harmonization.pan_ukbb.pan_ukbb_allele_frequency_panel_task import (
    PAN_UKBB_MANIFEST_READ_SPEC,
)

PAN_UKBB_VARIANT_MANIFEST = DownloadFileTask(
    meta=ReferenceFileMeta(
        id="pan_ukbb_variant_manifest",
        group="pan_ukbb",
        sub_group="variant_manifest",
        sub_folder=PurePath("raw"),
        extension=".txt.bgz",
        read_spec=PAN_UKBB_MANIFEST_READ_SPEC,
    ),
    url=(
        "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/"
        "full_variant_qc_metrics.txt.bgz"
    ),
    md5_hash="e70ebc8289f762dd8d5086f54e766654",
)
