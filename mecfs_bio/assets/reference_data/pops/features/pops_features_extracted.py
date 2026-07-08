"""
Task to gunzip the full POPs feature collection into a single feature table.

The resulting FileAsset is staged into a feature directory and munged by
PopsMungeFeatureDirectoryTask.
"""

from mecfs_bio.assets.reference_data.pops.features.pops_features_raw import (
    POPS_FEATURES_RAW,
)
from mecfs_bio.build_system.task.extract_gzip_task import ExtractGzipTextFileTask

POPS_FEATURES_EXTRACTED = ExtractGzipTextFileTask.create_for_reference_file(
    source_file_task=POPS_FEATURES_RAW,
    asset_id="pops_features_extracted",
)
