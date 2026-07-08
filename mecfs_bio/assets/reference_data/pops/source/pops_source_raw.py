"""
Task to download the POPs source code.

POPs (https://github.com/FinucaneLab/pops) is not distributed as a Python package
(no PyPI release, no setup.py), so we fetch a pinned GitHub source tarball and run
its scripts as subprocesses. The tarball also bundles the canonical gene-annotation
and control-feature files POPs uses, plus example data used by the system test.

The code is licensed GPLv3; we invoke it as a separate process (we do not vendor or
link it), so there is no license obligation on this repository's own code.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

# Pinned commit on FinucaneLab/pops master. The archive extracts to a single
# top-level directory named pops-<POPS_SOURCE_COMMIT>.
POPS_SOURCE_COMMIT = "76eb86cba10254490003c8f4dc7ff5ce492d3667"
POPS_SOURCE_TARBALL_TOP_DIR = f"pops-{POPS_SOURCE_COMMIT}"

POPS_SOURCE_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="pops",
        sub_group="source",
        sub_folder=PurePath("raw"),
        extension=".tar.gz",
        id="pops_source_raw",
    ),
    url=f"https://github.com/FinucaneLab/pops/archive/{POPS_SOURCE_COMMIT}.tar.gz",
    md5_hash="543f30363511b5de653a19dbe3dc0aa8",
)
