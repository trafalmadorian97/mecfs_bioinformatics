from pathlib import Path

from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.build_system.wf.synapse_downloader import SynapseDownloader
from mecfs_bio.util.file_io.ranged_tar import RangeReader


class _BytesRangeReader:
    def __init__(self, data: bytes) -> None:
        self.data = data

    def read_range(self, start: int, length: int) -> bytes:
        return self.data[start : start + length]


class _RecordingDownloader(SynapseDownloader):
    def __init__(self) -> None:
        self.calls: list[tuple[str, Path, str]] = []

    def download(self, synid: str, dest_dir: Path, expected_name: str) -> Path:
        self.calls.append((synid, dest_dir, expected_name))
        out = dest_dir / expected_name
        out.write_text("payload")
        return out

    def open_range_reader(self, synid: str) -> RangeReader:
        return _BytesRangeReader(b"header-bytes")


def test_wf_delegates_download_from_synapse(tmp_path: Path) -> None:
    fake = _RecordingDownloader()
    wf = make_wf(synapse_downloader=fake)

    out = wf.download_from_synapse("syn123", tmp_path, "protein.tar")

    assert out == tmp_path / "protein.tar"
    assert out.read_text() == "payload"
    assert fake.calls == [("syn123", tmp_path, "protein.tar")]


def test_wf_delegates_fetch_synapse_file_head() -> None:
    fake = _RecordingDownloader()
    wf = make_wf(synapse_downloader=fake)

    head = wf.fetch_synapse_file_head("syn123", 6)

    assert head == b"header"
