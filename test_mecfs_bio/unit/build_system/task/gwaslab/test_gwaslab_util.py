from pathlib import Path

from attrs import define

from mecfs_bio.build_system.task.gwaslab.gwaslab_util import (
    expected_reference_md5,
    gwaslab_download_ref_if_missing,
)

_REF = "1kg_eur_hg19"
_STALE_MD5 = "734069d895009d38c2f962bfbb6fab52"


@define
class FakeGwaslabRefStore:
    """Stands in for gwaslab's on-disk reference cache plus its config registry.

    gwaslab only registers a path once the file passes its own checksum test, and it
    skips the fetch entirely when a file of the expected name is already there unless
    asked to overwrite. This reproduces both behaviours.
    """

    expected_md5: str
    file_md5: str | None = None
    registered: Path | None = None
    downloads: list[bool] = []

    def path_lookup(self, ref: str) -> Path | None:
        return self.registered

    def download(self, ref: str, overwrite: bool) -> None:
        self.downloads.append(overwrite)
        if self.file_md5 is None or overwrite:
            self.file_md5 = self.expected_md5
        if self.file_md5 == self.expected_md5:
            self.registered = Path(f"/fake/.gwaslab/{ref}")

    def checksum(self, path: Path) -> str:
        assert self.file_md5 is not None
        return self.file_md5


def _make_store(file_md5: str | None = None) -> FakeGwaslabRefStore:
    """A store whose cached file has file_md5, or which has no cached file at all."""
    expected = expected_reference_md5(_REF)
    assert expected is not None
    return FakeGwaslabRefStore(expected_md5=expected, file_md5=file_md5, downloads=[])


def _good_file_store() -> FakeGwaslabRefStore:
    store = _make_store()
    store.file_md5 = store.expected_md5
    return store


def test_stale_registry_and_stale_file_forces_a_fresh_download():
    """The registry points nowhere (a config reset) while a superseded copy of the
    reference sits under the name gwaslab downloads to."""
    store = _make_store(file_md5=_STALE_MD5)

    result = gwaslab_download_ref_if_missing(
        _REF,
        path_lookup=store.path_lookup,
        downloader=store.download,
        checksum=store.checksum,
    )

    assert result == Path(f"/fake/.gwaslab/{_REF}")
    assert store.downloads == [False, True]


def test_stale_registry_with_a_good_file_registers_without_overwriting():
    store = _good_file_store()

    result = gwaslab_download_ref_if_missing(
        _REF,
        path_lookup=store.path_lookup,
        downloader=store.download,
        checksum=store.checksum,
    )

    assert result == Path(f"/fake/.gwaslab/{_REF}")
    assert store.downloads == [False]


def test_registered_reference_with_a_matching_checksum_is_not_downloaded():
    store = _good_file_store()
    store.registered = Path("/fake/.gwaslab/already_here")

    result = gwaslab_download_ref_if_missing(
        _REF,
        path_lookup=store.path_lookup,
        downloader=store.download,
        checksum=store.checksum,
    )

    assert result == Path("/fake/.gwaslab/already_here")
    assert store.downloads == []


def test_registered_reference_with_a_stale_checksum_is_redownloaded():
    store = _make_store(file_md5=_STALE_MD5)
    store.registered = Path("/fake/.gwaslab/already_here")

    result = gwaslab_download_ref_if_missing(
        _REF,
        path_lookup=store.path_lookup,
        downloader=store.download,
        checksum=store.checksum,
    )

    assert result == Path(f"/fake/.gwaslab/{_REF}")
    assert store.downloads == [False, True]
