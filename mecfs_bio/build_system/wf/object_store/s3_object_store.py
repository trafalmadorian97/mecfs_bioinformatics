"""Production ObjectStore backed by S3.

Not unit-tested: boto3 client calls here are only exercised through the
FakeObjectStore substitute in unit tests, plus a one-time manual maintainer
staging run against real S3 (see design spec). The exact shape of
get_object_attributes's Checksum payload is inferred from the boto3/S3 API
docs rather than verified against a live bucket, so treat the sha256 return
of head() as a best-effort convenience, not a guaranteed contract.

Composite-checksum warning: upload_fileobj uses boto3's TransferManager,
which multipart-uploads any file over roughly 8 MiB, and our reference files
(hundreds of GB) always take that path. For a multipart upload, the
ChecksumSHA256 that S3 stores and that get_object_attributes reports back is
a COMPOSITE checksum (a hash of the per-part checksums), not a plain
byte-for-byte SHA-256 of the object. Consequences for callers:

- An S3-reported sha256 (from head() or upload_from_url()) may only be
  compared against another S3-reported sha256 for the SAME object (e.g. a
  dedup check comparing head().sha256 against a value previously returned by
  upload_from_url()). Both sides go through the same composite scheme, so
  the comparison is meaningful.
- An S3-reported sha256 must NEVER be compared against an independently
  computed local digest (e.g. `sha256sum` on the source file before upload);
  for a multipart object those two values are not the same hash and a
  mismatch does not indicate corruption, nor does a match guarantee byte
  equality.
- Forcing whole-object checksums (S3's ChecksumType=FULL_OBJECT) would avoid
  this, but is deliberately not done here; whether it is worth the extra
  configuration is left for validation during the real staging run, not
  decided speculatively in this untested code path.
"""

import math
from collections.abc import Callable
from typing import Any
from urllib.request import urlopen

import boto3
from attrs import Factory, frozen
from boto3.s3.transfer import TransferConfig
from botocore.exceptions import ClientError

from mecfs_bio.build_system.wf.object_store.base_object_store import (
    ObjectHead,
    ObjectStore,
)

_NOT_FOUND_ERROR_CODES = {"404", "NoSuchKey", "NotFound"}
_URLOPEN_TIMEOUT_SECONDS = 60

# s3transfer's default multipart part size. A multipart upload may have at most
# _S3_MAX_PARTS parts, so this default caps an object at 10,000 * 8 MiB = 80 GiB;
# a larger object overflows the part count and S3 rejects the UploadPart call.
_DEFAULT_MULTIPART_CHUNKSIZE = 8 * 1024**2
# S3's hard limit on parts per multipart upload.
_S3_MAX_PARTS = 10_000
# Aim well under the hard limit so a slightly-larger-than-reported source, or any
# off-by-one in s3transfer's own part accounting, cannot push us over.
_MULTIPART_PART_BUDGET = 8_000
_MIB = 1024**2


def multipart_chunksize_for_size(total_bytes: int) -> int:
    """Choose an S3 multipart chunk size that fits total_bytes within the part limit.

    S3 allows at most _S3_MAX_PARTS parts per multipart upload, so the default
    _DEFAULT_MULTIPART_CHUNKSIZE only works up to about 80 GiB; a larger object needs
    a bigger part size or S3 rejects it. This returns the larger of the default and
    the smallest MiB-aligned chunk that splits total_bytes into at most
    _MULTIPART_PART_BUDGET parts, so objects at or below the default's capacity are
    left untouched while huge objects (the hundreds-of-GB LD matrices) stay well
    within the part limit.
    """
    if total_bytes <= 0:
        return _DEFAULT_MULTIPART_CHUNKSIZE
    needed = math.ceil(total_bytes / _MULTIPART_PART_BUDGET)
    mib_aligned = math.ceil(needed / _MIB) * _MIB
    return max(_DEFAULT_MULTIPART_CHUNKSIZE, mib_aligned)


def _split_s3_uri(uri: str) -> tuple[str, str]:
    """Split an s3://bucket/key uri into (bucket, key)."""
    assert uri.startswith("s3://"), f"Expected an s3:// uri, got {uri!r}"
    without_scheme = uri[len("s3://") :]
    bucket, _, key = without_scheme.partition("/")
    assert bucket and key, f"Expected s3://bucket/key, got {uri!r}"
    return bucket, key


@frozen
class S3ObjectStore(ObjectStore):
    """ObjectStore backed by a real S3 bucket via boto3.




    Uploads use ONEZONE_IA storage (this is scratch/staging data, not
    durable-tier) and request an S3-managed SHA-256 checksum so head() can
    later report it without re-downloading the object. See the module
    docstring for the composite-checksum caveat on large (multipart)
    objects: the returned sha256 is only safe to compare against another
    S3-reported sha256, never against an independently computed local
    digest.
    """

    # Typed Any because boto3's S3 client class is generated at runtime; the only
    # precise static type is mypy_boto3_s3.S3Client from boto3-stubs, not a dependency
    # here. The name has no leading underscore: attrs' underscore-stripping for init
    # args has tripped up the ty typechecker.
    client: Any = Factory(lambda: boto3.client("s3"))

    def head(self, uri: str) -> ObjectHead | None:
        """Return size and S3-reported checksum for uri, or None if absent.

        The returned sha256 (if not None) is S3's checksum for the stored
        object; for a multipart-uploaded object (see module docstring) that
        is a composite checksum, not a plain SHA-256 of the bytes. Only
        compare it against another S3-reported checksum.
        """
        bucket, key = _split_s3_uri(uri)
        try:
            head_response = self.client.head_object(Bucket=bucket, Key=key)
        except ClientError as error:
            code = error.response.get("Error", {}).get("Code")
            if code in _NOT_FOUND_ERROR_CODES:
                return None
            raise
        size_bytes = head_response["ContentLength"]

        attributes_response = self.client.get_object_attributes(
            Bucket=bucket, Key=key, ObjectAttributes=["Checksum"]
        )
        checksum = attributes_response.get("Checksum") or {}
        sha256 = checksum.get("ChecksumSHA256")

        return ObjectHead(size_bytes=size_bytes, sha256=sha256)

    def upload_from_url(
        self,
        source_url: str,
        uri: str,
        on_progress: Callable[[int], None] | None = None,
    ) -> str:
        """Stream source_url's body into uri and return S3's stored sha256.
        Summary:
          urlopen returns an http.client.HTTPResponse, which is a file-like object: it exposes .read(n) that pulls the
         next n bytes off the underlying TCP socket on demand. upload_fileobj pumps that file-like-object through a bounded,
        concurrent S3 multipart upload

        The returned string is S3's checksum for the stored object; for a
        multipart upload (see module docstring, and note that our reference
        files are always multipart) that is a composite checksum, not a
        plain SHA-256 of source_url's bytes. Only compare it against another
        S3-reported checksum (e.g. a later head() call), never against an
        independently computed local digest of the source.

        on_progress is forwarded to boto3 as the upload Callback: s3transfer
        invokes it from its worker threads with the bytes moved per part, so a
        supplied callback must be thread-safe. None (the default) means boto3
        reports no progress.
        """
        bucket, key = _split_s3_uri(uri)
        with urlopen(source_url, timeout=_URLOPEN_TIMEOUT_SECONDS) as response:  # noqa: S310
            # response.length is the remaining Content-Length before any read, i.e.
            # the object size (None for a chunked source, which falls back to the
            # default). Size the multipart chunk from it so a hundreds-of-GB object
            # does not overflow S3's 10,000-part limit at the default 8 MiB part size.
            chunksize = multipart_chunksize_for_size(response.length or 0)
            self.client.upload_fileobj(
                response,
                bucket,
                key,
                ExtraArgs={
                    "StorageClass": "ONEZONE_IA",
                    "ChecksumAlgorithm": "SHA256",
                },
                Config=TransferConfig(multipart_chunksize=chunksize),
                Callback=on_progress,
            )

        attributes_response = self.client.get_object_attributes(
            Bucket=bucket, Key=key, ObjectAttributes=["Checksum"]
        )
        checksum = attributes_response.get("Checksum") or {}
        sha256 = checksum.get("ChecksumSHA256")
        assert sha256, (
            f"Expected S3 to report a stored SHA-256 checksum for {uri} after "
            f"upload_fileobj with ChecksumAlgorithm=SHA256"
        )
        return sha256
