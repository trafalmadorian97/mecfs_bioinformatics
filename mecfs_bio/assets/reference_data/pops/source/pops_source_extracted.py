"""
Task to extract the POPs source tarball to a directory.

Strips the pops-<commit> top-level directory the GitHub archive wraps everything
in, so the resulting DirectoryAsset root directly contains pops.py,
munge_feature_directory.py, and the example directory.
"""

from mecfs_bio.assets.reference_data.pops.source.pops_source_raw import (
    POPS_SOURCE_RAW,
    POPS_SOURCE_TARBALL_TOP_DIR,
)
from mecfs_bio.build_system.task.extract_tar_gzip_task import ExtractTarGzipTask

POPS_SOURCE_EXTRACTED = ExtractTarGzipTask.create(
    asset_id="pops_source_extracted",
    source_task=POPS_SOURCE_RAW,
    sub_folder_name_inside_tar=POPS_SOURCE_TARBALL_TOP_DIR,
)
