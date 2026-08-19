import math

from mecfs_bio.build_system.wf.object_store.s3_object_store import (
    _DEFAULT_MULTIPART_CHUNKSIZE,
    multipart_chunksize_for_size,
)

# S3's hard limit on the number of parts in a single multipart upload.
_S3_HARD_PART_LIMIT = 10_000


def test_huge_object_stays_within_s3_part_limit() -> None:
    # The 192 GiB LD matrix that failed in production: at the 8 MiB default chunk it
    # needs ~24,600 parts and S3 rejects part 10,001. The chosen chunk must keep the
    # part count within S3's limit.
    total_bytes = int(192.38 * 1024**3)
    chunksize = multipart_chunksize_for_size(total_bytes)
    parts = math.ceil(total_bytes / chunksize)
    assert parts <= _S3_HARD_PART_LIMIT


def test_object_within_default_capacity_keeps_the_default_chunk() -> None:
    # A file the 8 MiB default already handles (well under 80 GiB) must not be
    # perturbed, so ordinary uploads keep boto3's default behavior.
    ten_gib = 10 * 1024**3
    assert multipart_chunksize_for_size(ten_gib) == _DEFAULT_MULTIPART_CHUNKSIZE
