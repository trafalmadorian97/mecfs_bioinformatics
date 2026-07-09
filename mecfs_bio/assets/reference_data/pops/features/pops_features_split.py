"""
Task to split the monolithic POPs feature collection into a directory of column
chunks.

This is the memory-safe producer that feeds PopsMungeFeatureDirectoryTask: the full
PoPS.features file is a single ~21.5GB table that the munge step cannot load whole,
so it is streamed into manageable per-chunk files first. See
PopsSplitFeatureFileTask.
"""

from mecfs_bio.assets.reference_data.pops.features.pops_features_extracted import (
    POPS_FEATURES_EXTRACTED,
)
from mecfs_bio.build_system.task.pops.pops_split_features_task import (
    PopsSplitFeatureFileTask,
)

POPS_FEATURES_SPLIT = PopsSplitFeatureFileTask.create(
    asset_id="pops_features_split",
    source_features_task=POPS_FEATURES_EXTRACTED,
)
