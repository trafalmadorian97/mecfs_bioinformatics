"""
Generic access to objects in S3-compatible buckets, for Tasks whose inputs live there.

A Task holds an S3ObjectRef (endpoint, bucket, key, and the name of the AWS credentials
profile to sign with), never keys or presigned URLs. The WF's S3Access turns a ref into object
metadata (HEAD) or a time-limited presigned GET URL, which the WF's downloaders then fetch.

Credentials come from the named profile in the AWS shared credentials file. If the keys
behind a profile are rotated or expire mid-run, the first 403 drops the cached client and
rebuilds it (re-reading the profile); a second 403 fails with a message saying how to fix the
profile.
"""

import base64
import threading
from abc import ABC, abstractmethod
from collections.abc import Callable

import boto3
from attrs import define, field, frozen
from botocore.config import Config
from botocore.exceptions import ClientError, ProfileNotFound

S3_CREDENTIALS_HELP = (
    "S3 access with AWS profile {profile!r} was refused or the profile is missing. "
    "Set its keys with: aws configure --profile {profile}"
)
_FORBIDDEN = 403
_MD5_METADATA_KEY = "md5chksum"


@frozen(slots=True)
class S3Endpoint:
    """
    url: the S3 API endpoint, e.g. https://s3-ext.decode.is.
    region: the signing region.
    credentials_profile: the AWS shared-credentials profile whose keys sign requests.
    """

    url: str
    region: str
    credentials_profile: str


@frozen(slots=True)
class S3ObjectRef:
    endpoint: S3Endpoint
    bucket: str
    key: str


@frozen(slots=True)
class S3ObjectInfo:
    """
    size: object size in bytes.
    md5: lowercase hex md5 from the base64 x-amz-meta-md5chksum metadata, or None if absent.
    etag: the ETag without its quotes (the md5 for objects not uploaded in parts).
    """

    size: int
    md5: str | None
    etag: str


class S3CredentialsError(Exception):
    pass


class S3Client(ABC):
    """The few S3 operations S3Access needs, so tests can supply a fake."""

    @abstractmethod
    def head_object(self, bucket: str, key: str) -> dict:
        pass

    @abstractmethod
    def generate_presigned_url(self, bucket: str, key: str, expires_s: int) -> str:
        pass


class S3ClientFactory(ABC):
    @abstractmethod
    def make(self, endpoint: S3Endpoint) -> S3Client:
        pass


@frozen(slots=True)
class _Boto3Client(S3Client):
    client: object

    def head_object(self, bucket: str, key: str) -> dict:
        return self.client.head_object(Bucket=bucket, Key=key)  # ty: ignore[missing-attribute]

    def generate_presigned_url(self, bucket: str, key: str, expires_s: int) -> str:
        return self.client.generate_presigned_url(  # ty: ignore[missing-attribute]
            "get_object", Params={"Bucket": bucket, "Key": key}, ExpiresIn=expires_s
        )


@frozen(slots=True)
class Boto3ClientFactory(S3ClientFactory):
    def make(self, endpoint: S3Endpoint) -> S3Client:
        session = boto3.Session(profile_name=endpoint.credentials_profile)
        return _Boto3Client(
            session.client(
                "s3",
                endpoint_url=endpoint.url,
                region_name=endpoint.region,
                config=Config(
                    signature_version="s3v4", s3={"addressing_style": "path"}
                ),
            )
        )


class S3Access(ABC):
    @abstractmethod
    def head(self, ref: S3ObjectRef) -> S3ObjectInfo:
        pass

    @abstractmethod
    def presigned_url(self, ref: S3ObjectRef, expires_s: int) -> str:
        pass


def _is_forbidden(error: ClientError) -> bool:
    status = error.response.get("ResponseMetadata", {}).get("HTTPStatusCode")
    return status == _FORBIDDEN or error.response.get("Error", {}).get("Code") in (
        "403",
        "Forbidden",
        "AccessDenied",
    )


def _object_info(head: dict) -> S3ObjectInfo:
    encoded = head.get("Metadata", {}).get(_MD5_METADATA_KEY)
    return S3ObjectInfo(
        size=int(head["ContentLength"]),
        md5=None if encoded is None else base64.b64decode(encoded).hex(),
        etag=str(head["ETag"]).strip('"'),
    )


@define(slots=True)
class Boto3S3Access(S3Access):
    """S3Access through one cached client per endpoint, rebuilt once after a 403."""

    client_factory: S3ClientFactory = field(factory=Boto3ClientFactory)
    _clients: dict[S3Endpoint, S3Client] = field(factory=dict, init=False)
    _lock: threading.Lock = field(factory=threading.Lock, init=False)

    def head(self, ref: S3ObjectRef) -> S3ObjectInfo:
        return _object_info(
            self._with_client(
                ref.endpoint, lambda c: c.head_object(ref.bucket, ref.key)
            )
        )

    def presigned_url(self, ref: S3ObjectRef, expires_s: int) -> str:
        return self._with_client(
            ref.endpoint,
            lambda c: c.generate_presigned_url(ref.bucket, ref.key, expires_s),
        )

    def _client(self, endpoint: S3Endpoint, rebuild: bool) -> S3Client:
        with self._lock:
            if rebuild or endpoint not in self._clients:
                try:
                    self._clients[endpoint] = self.client_factory.make(endpoint)
                except ProfileNotFound as e:
                    raise S3CredentialsError(_help(endpoint)) from e
            return self._clients[endpoint]

    def _with_client[T](self, endpoint: S3Endpoint, op: Callable[[S3Client], T]) -> T:
        try:
            return op(self._client(endpoint, rebuild=False))
        except ClientError as e:
            if not _is_forbidden(e):
                raise
        try:
            return op(self._client(endpoint, rebuild=True))
        except ClientError as e:
            if _is_forbidden(e):
                raise S3CredentialsError(_help(endpoint)) from e
            raise


def _help(endpoint: S3Endpoint) -> str:
    return S3_CREDENTIALS_HELP.format(profile=endpoint.credentials_profile)
