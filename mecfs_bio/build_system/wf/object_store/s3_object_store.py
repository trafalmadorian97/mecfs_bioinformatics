"""Production ObjectStore backed by S3.

Not unit-tested: boto3 client calls here are only exercised through the
FakeObjectStore substitute in unit tests, plus a one-time manual maintainer
staging run against real S3 (see design spec). The exact shape of
get_object_attributes's Checksum payload is inferred from the boto3/S3 API
docs rather than verified against a live bucket, so treat the sha256 return
of head() as a best-effort convenience, not a guaranteed contract.
"""

from typing import Any
from urllib.request import urlopen

import boto3
from botocore.exceptions import ClientError

from mecfs_bio.build_system.wf.object_store.base_object_store import (
    ObjectHead,
    ObjectStore,
)

_NOT_FOUND_ERROR_CODES = {"404", "NoSuchKey", "NotFound"}


def _split_s3_uri(uri: str) -> tuple[str, str]:
    """Split an s3://bucket/key uri into (bucket, key)."""
    assert uri.startswith("s3://"), f"Expected an s3:// uri, got {uri!r}"
    without_scheme = uri[len("s3://") :]
    bucket, _, key = without_scheme.partition("/")
    assert bucket and key, f"Expected s3://bucket/key, got {uri!r}"
    return bucket, key


class S3ObjectStore(ObjectStore):
    """ObjectStore backed by a real S3 bucket via boto3.

    Uploads use ONEZONE_IA storage (this is scratch/staging data, not
    durable-tier) and request an S3-managed SHA-256 checksum so head() can
    later report it without re-downloading the object.
    """

    def __init__(self, client: Any = None) -> None:
        self._client = client if client is not None else boto3.client("s3")

    def head(self, uri: str) -> ObjectHead | None:
        bucket, key = _split_s3_uri(uri)
        try:
            head_response = self._client.head_object(Bucket=bucket, Key=key)
        except ClientError as error:
            code = error.response.get("Error", {}).get("Code")
            if code in _NOT_FOUND_ERROR_CODES:
                return None
            raise
        size_bytes = head_response["ContentLength"]

        attributes_response = self._client.get_object_attributes(
            Bucket=bucket, Key=key, ObjectAttributes=["Checksum"]
        )
        checksum = attributes_response.get("Checksum") or {}
        sha256 = checksum.get("ChecksumSHA256")

        return ObjectHead(size_bytes=size_bytes, sha256=sha256)

    def upload_from_url(self, source_url: str, uri: str) -> str:
        bucket, key = _split_s3_uri(uri)
        with urlopen(source_url) as response:  # noqa: S310
            self._client.upload_fileobj(
                response,
                bucket,
                key,
                ExtraArgs={
                    "StorageClass": "ONEZONE_IA",
                    "ChecksumAlgorithm": "SHA256",
                },
            )

        attributes_response = self._client.get_object_attributes(
            Bucket=bucket, Key=key, ObjectAttributes=["Checksum"]
        )
        checksum = attributes_response.get("Checksum") or {}
        sha256 = checksum.get("ChecksumSHA256")
        assert sha256, (
            f"Expected S3 to report a stored SHA-256 checksum for {uri} after "
            f"upload_fileobj with ChecksumAlgorithm=SHA256"
        )
        return sha256
