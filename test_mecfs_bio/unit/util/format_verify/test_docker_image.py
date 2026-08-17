import pytest

from mecfs_bio.util.format_verify.docker_image import is_valid_docker_image_reference

_VALID = [
    "ghcr.io/trafalmadorian97/gctb:2.5.5",
    "debian:stable-slim",
    "zhiliz/sbayesrc",
    "ubuntu",
    "busybox:latest",
    "localhost:5000/foo",
    "gctb:test",
    "ghcr.io/org/repo@sha256:" + "a" * 64,
    "repo@sha256:" + "0" * 64,
]

_INVALID = [
    "",
    "   ",
    "has space",
    "gctb:",
    "Upper/Case",
    "UPPERCASE",
    "bad_/name",
    "repo@sha256:xyz",
    "repo@notadigest",
    "/leadingslash",
]


@pytest.mark.parametrize("ref", _VALID)
def test_accepts_standard_references(ref: str) -> None:
    assert is_valid_docker_image_reference(ref)


@pytest.mark.parametrize("ref", _INVALID)
def test_rejects_malformed_references(ref: str) -> None:
    assert not is_valid_docker_image_reference(ref)
