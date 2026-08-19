import pytest

from mecfs_bio.util.format_verify.s3_uri import (
    is_valid_s3_bucket_name,
    is_valid_s3_uri,
)

_VALID_URIS = [
    "s3://mecfs-bio-remote-exec-scratch/remote-exec-scratch",
    "s3://bucket",
    "s3://bucket/",
    "s3://bucket/pre/fix/",
    "s3://my.dotted.bucket/key",
    "s3://mecfs-bio-reference-data/sbayesrc/reference/Imputed13M/v1/",
]

_INVALID_URIS = [
    "",
    "   ",
    "s3://bucket/has space",
    "mecfs-bio-reference-data/no-scheme",
    "https://bucket/key",
    "s3://Bad_Bucket/key",
    "s3://bad..bucket/key",
    "s3://ab/key",
    "s3://192.168.0.1/key",
    "s3:///key",
]


@pytest.mark.parametrize("uri", _VALID_URIS)
def test_accepts_valid_uris(uri: str) -> None:
    assert is_valid_s3_uri(uri)


@pytest.mark.parametrize("uri", _INVALID_URIS)
def test_rejects_invalid_uris(uri: str) -> None:
    assert not is_valid_s3_uri(uri)


@pytest.mark.parametrize(
    "name",
    ["bucket", "abc", "my.dotted.bucket", "mecfs-bio-remote-exec-scratch"],
)
def test_accepts_valid_bucket_names(name: str) -> None:
    assert is_valid_s3_bucket_name(name)


@pytest.mark.parametrize(
    "name",
    ["ab", "a" * 64, "Bad_Bucket", "bad..bucket", "192.168.0.1", "-lead", "trail-"],
)
def test_rejects_invalid_bucket_names(name: str) -> None:
    assert not is_valid_s3_bucket_name(name)
