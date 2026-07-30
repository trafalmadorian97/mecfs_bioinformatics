from pathlib import Path

from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.build_system.wf.synapse_downloader import SynapseDownloader


class _RecordingDownloader(SynapseDownloader):
    def __init__(self) -> None:
        self.calls: list[tuple[str, Path, str]] = []
        self.head_calls: list[tuple[str, int]] = []

    def download(self, synid: str, dest_dir: Path, expected_name: str) -> Path:
        self.calls.append((synid, dest_dir, expected_name))
        out = dest_dir / expected_name
        out.write_text("payload")
        return out

    def fetch_file_head(self, synid: str, n_bytes: int) -> bytes:
        self.head_calls.append((synid, n_bytes))
        return b"header-bytes"[:n_bytes]


def test_wf_delegates_download_from_synapse(tmp_path: Path):
    fake = _RecordingDownloader()
    wf = make_wf(synapse_downloader=fake)

    out = wf.download_from_synapse("syn123", tmp_path, "protein.tar")

    assert out == tmp_path / "protein.tar"
    assert out.read_text() == "payload"
    assert fake.calls == [("syn123", tmp_path, "protein.tar")]


def test_wf_delegates_fetch_synapse_file_head():
    fake = _RecordingDownloader()
    wf = make_wf(synapse_downloader=fake)

    head = wf.fetch_synapse_file_head("syn123", 6)

    assert head == b"header"
    assert fake.head_calls == [("syn123", 6)]
