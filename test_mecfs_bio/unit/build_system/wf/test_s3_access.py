import base64

import pytest
from attrs import define, field
from botocore.exceptions import ClientError

from mecfs_bio.build_system.wf.s3_access import (
    Boto3S3Access,
    S3Client,
    S3ClientFactory,
    S3CredentialsError,
    S3Endpoint,
    S3ObjectRef,
)

_ENDPOINT = S3Endpoint(
    url="https://s3.example", region="us-east-1", credentials_profile="p"
)
_OTHER_ENDPOINT = S3Endpoint(
    url="https://other.example", region="us-east-1", credentials_profile="p"
)
_MD5_BYTES = bytes.fromhex("933a" + "00" * 13 + "ff")


def _forbidden() -> ClientError:
    return ClientError(
        {
            "Error": {"Code": "403", "Message": "Forbidden"},
            "ResponseMetadata": {"HTTPStatusCode": 403},
        },
        "HeadObject",
    )


@define(slots=True)
class _FakeClient(S3Client):
    forbidden: bool = False
    metadata: dict[str, str] = field(factory=dict)

    def head_object(self, bucket: str, key: str) -> dict:
        if self.forbidden:
            raise _forbidden()
        return {"ContentLength": 10, "ETag": '"abc"', "Metadata": self.metadata}

    def generate_presigned_url(self, bucket: str, key: str, expires_s: int) -> str:
        return f"https://signed/{bucket}/{key}"


@define(slots=True)
class _ScriptedFactory(S3ClientFactory):
    clients: list[_FakeClient]
    made: int = 0

    def make(self, endpoint: S3Endpoint) -> S3Client:
        client = self.clients[min(self.made, len(self.clients) - 1)]
        self.made += 1
        return client


def _ref(endpoint: S3Endpoint = _ENDPOINT) -> S3ObjectRef:
    return S3ObjectRef(endpoint=endpoint, bucket="b", key="k.txt.gz")


def test_client_reused_per_endpoint() -> None:
    factory = _ScriptedFactory(clients=[_FakeClient()])
    access = Boto3S3Access(client_factory=factory)
    access.head(_ref())
    access.head(_ref())
    assert factory.made == 1
    access.head(_ref(_OTHER_ENDPOINT))
    assert factory.made == 2


def test_403_rebuilds_client_once_then_succeeds() -> None:
    factory = _ScriptedFactory(clients=[_FakeClient(forbidden=True), _FakeClient()])
    access = Boto3S3Access(client_factory=factory)
    info = access.head(_ref())
    assert info.size == 10
    assert factory.made == 2


def test_persistent_403_raises_credentials_error() -> None:
    factory = _ScriptedFactory(clients=[_FakeClient(forbidden=True)])
    access = Boto3S3Access(client_factory=factory)
    with pytest.raises(S3CredentialsError):
        access.head(_ref())


def test_head_decodes_md5_metadata() -> None:
    encoded = base64.b64encode(_MD5_BYTES).decode()
    with_md5 = Boto3S3Access(
        client_factory=_ScriptedFactory(
            clients=[_FakeClient(metadata={"md5chksum": encoded})]
        )
    ).head(_ref())
    assert with_md5.md5 == _MD5_BYTES.hex()

    without_md5 = Boto3S3Access(
        client_factory=_ScriptedFactory(clients=[_FakeClient()])
    ).head(_ref())
    assert without_md5.md5 is None
    assert without_md5.etag == "abc"
