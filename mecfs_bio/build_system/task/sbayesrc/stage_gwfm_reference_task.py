"""
Stage the pinned GWFM LD reference bundle into an S3 bucket.

This Task uploads each file in GWFM_REFERENCE_BUNDLE (Task 1) to the bucket, skipping
files that are already present with a matching size (and checksum, once pinned), and
writes a small JSON marker asset recording that staging happened.

The marker content is derived only from the pinned bundle and the target bucket, never
from an S3-reported checksum. S3-reported checksums are, for multipart uploads (the
206GB LD matrix file), composite checksums of per-part hashes rather than a byte-for-byte
SHA-256, and they can differ across machines and across re-uploads of the same bytes. If
the marker embedded a live S3 checksum, the build system's trace (which hashes asset
content) would be non-deterministic across collaborators and re-runs, forcing endless
rebuilds. Dedup therefore only ever compares one S3-reported value (a stored object's
head) against another S3-reported value (the pinned constant, once recorded from a prior
run's logs).
"""

import json
import threading
from collections.abc import Callable
from pathlib import Path, PurePath

import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.sbayesrc.gctb_gwfm_constants import (
    GWFM_REFERENCE_BUNDLE,
    GWFM_REFERENCE_VERSION,
    MARKER_FILES_KEY,
    MARKER_PREFIX_KEY,
    MARKER_VERSION_KEY,
    gwfm_reference_prefix,
)
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()

_MARKER_ASSET_ID = "gwfm_reference_marker"

# Log upload progress at most once per this fraction of the file, so the ~192 GiB LD
# file yields a couple dozen lines rather than one per 8 MiB part.
_PROGRESS_LOG_STEP_FRACTION = 0.05


def make_upload_progress_logger(
    filename: str,
    total_bytes: int,
    log: Callable[..., object] = logger.info,
) -> Callable[[int], None]:
    """Build a thread-safe on_progress callback that logs throttled upload percentage.

    The object store invokes the returned callback from several worker threads, each
    time with the number of bytes just transferred, so it accumulates under a lock and
    emits a line only when the completed fraction crosses the next
    _PROGRESS_LOG_STEP_FRACTION boundary (and once at completion). total_bytes is the
    file's known size, used to turn the byte deltas into a percentage. log is injected
    so tests can capture the calls.
    """
    lock = threading.Lock()
    seen = 0
    next_threshold = _PROGRESS_LOG_STEP_FRACTION

    def on_progress(bytes_transferred: int) -> None:
        nonlocal seen, next_threshold
        with lock:
            seen += bytes_transferred
            fraction = seen / total_bytes if total_bytes else 1.0
            if fraction < next_threshold and seen < total_bytes:
                return
            log(
                "staging progress",
                filename=filename,
                percent=round(fraction * 100),
                transferred_gib=round(seen / 1024**3, 2),
                total_gib=round(total_bytes / 1024**3, 2),
            )
            # Advance past every band this delta crossed so a single large delta logs
            # once, not once per band.
            while next_threshold <= fraction:
                next_threshold += _PROGRESS_LOG_STEP_FRACTION

    return on_progress


@frozen
class StageGwfmReferenceTask(Task):
    """
    Stages the pinned GWFM reference bundle into an S3 bucket, skipping files
    already present with a matching size and checksum, and produces a deterministic
    JSON marker asset recording the staged S3 prefix and file list.
    """

    meta: Meta
    bucket: str

    @property
    def deps(self) -> list["Task"]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> FileAsset:
        assert GWFM_REFERENCE_BUNDLE, "GWFM reference bundle is empty"

        ref_prefix = gwfm_reference_prefix(GWFM_REFERENCE_VERSION)
        for bundle_file in GWFM_REFERENCE_BUNDLE:
            uri = f"s3://{self.bucket}/{ref_prefix}{bundle_file.filename}"
            head = wf.object_store.head(uri)
            needs_upload = (
                head is None
                or head.size_bytes != bundle_file.size_bytes
                or (
                    bundle_file.sha256 is not None and head.sha256 != bundle_file.sha256
                )
            )
            if needs_upload:
                # The returned sha256 is S3-reported (possibly a composite, multipart
                # checksum); it is only ever compared against another S3-reported
                # value, never against a locally computed digest.
                logger.info(
                    f"Uploading GFWM reference file, filename: {bundle_file.filename}, source url: {bundle_file.source_url}, uri :{uri}"
                )
                sha256 = wf.object_store.upload_from_url(
                    source_url=bundle_file.source_url,
                    uri=uri,
                    on_progress=make_upload_progress_logger(
                        filename=bundle_file.filename,
                        total_bytes=bundle_file.size_bytes,
                    ),
                )
                logger.info(
                    "uploaded GWFM reference file",
                    filename=bundle_file.filename,
                    uri=uri,
                    sha256=sha256,
                )
            else:
                logger.info(
                    "GWFM reference file already staged",
                    filename=bundle_file.filename,
                    uri=uri,
                    sha256=head.sha256,
                )

        marker = {
            MARKER_VERSION_KEY: GWFM_REFERENCE_VERSION,
            MARKER_PREFIX_KEY: f"s3://{self.bucket}/{ref_prefix}",
            MARKER_FILES_KEY: sorted(f.filename for f in GWFM_REFERENCE_BUNDLE),
        }
        target = scratch_dir / f"{_MARKER_ASSET_ID}.json"
        target.write_text(json.dumps(marker, sort_keys=True))
        return FileAsset(target)

    @classmethod
    def create(cls, bucket: str) -> "StageGwfmReferenceTask":
        meta = ReferenceFileMeta(
            group="sbayesrc_gwfm",
            sub_group="ld_reference",
            sub_folder=PurePath(GWFM_REFERENCE_VERSION),
            id=AssetId(_MARKER_ASSET_ID),
            extension=".json",
        )
        return cls(meta=meta, bucket=bucket)
