"""
Task to munge the full POPs feature collection into the chunked matrix format
pops.py consumes.

This is a single shared asset: the munge step depends only on the POPs source and
the raw feature collection, not on any GWAS, so the (multi-gigabyte, slow) munge is
run once here and reused across every per-GWAS POPs run rather than re-munged per
trait.
"""

from mecfs_bio.assets.reference_data.pops.features.pops_features_extracted import (
    POPS_FEATURES_EXTRACTED,
)
from mecfs_bio.assets.reference_data.pops.source.pops_source_extracted import (
    POPS_SOURCE_EXTRACTED,
)
from mecfs_bio.build_system.task.pops.pops_munge_task import (
    PopsMungeFeatureDirectoryTask,
)

POPS_FEATURES_MUNGED = PopsMungeFeatureDirectoryTask.create(
    asset_id="pops_features_munged",
    pops_source_task=POPS_SOURCE_EXTRACTED,
    raw_features_task=POPS_FEATURES_EXTRACTED,
)
