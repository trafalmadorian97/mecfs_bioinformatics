import io
import tarfile

import pytest

from mecfs_bio.util.file_io.ranged_tar import list_tar_members, read_tar_member

_SHORT_NAME = "protein/discovery_chr1_x.gz"
_LONG_NAME = "protein_" + "a" * 120 + "/discovery_chr2_" + "b" * 60 + ".gz"


class _BytesRangeReader:
    def __init__(self, data: bytes) -> None:
        self.data = data
        self.n_reads = 0

    def read_range(self, start: int, length: int) -> bytes:
        self.n_reads += 1
        return self.data[start : start + length]


def _tar_bytes(files: dict[str, bytes], tar_format: int) -> bytes:
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w", format=tar_format) as tar:
        directory = tarfile.TarInfo("protein")
        directory.type = tarfile.DIRTYPE
        tar.addfile(directory)
        for name, payload in files.items():
            info = tarfile.TarInfo(name)
            info.size = len(payload)
            tar.addfile(info, io.BytesIO(payload))
    return buffer.getvalue()


@pytest.mark.parametrize("tar_format", [tarfile.GNU_FORMAT, tarfile.USTAR_FORMAT])
def test_lists_regular_members_and_reads_their_exact_data(tar_format: int) -> None:
    # USTAR cannot hold the long name, so only GNU exercises long-name entries.
    files = {_SHORT_NAME: b"first payload", "protein/chr3.gz": b"x" * 5000}
    if tar_format == tarfile.GNU_FORMAT:
        files[_LONG_NAME] = b"second payload"
    reader = _BytesRangeReader(_tar_bytes(files, tar_format))

    members = list_tar_members(reader)

    assert [m.name for m in members] == list(files)
    assert {m.name: read_tar_member(reader, m) for m in members} == files


def test_listing_skips_member_data() -> None:
    files = {f"protein/chr{i}.gz": b"z" * 100_000 for i in range(5)}
    reader = _BytesRangeReader(_tar_bytes(files, tarfile.GNU_FORMAT))

    list_tar_members(reader)

    # One read per member header plus the closing zero block, never one per data block.
    assert reader.n_reads <= len(files) + 2


def test_pax_headers_fail_fast() -> None:
    reader = _BytesRangeReader(_tar_bytes({_LONG_NAME: b"p"}, tarfile.PAX_FORMAT))
    with pytest.raises(AssertionError):
        list_tar_members(reader)
