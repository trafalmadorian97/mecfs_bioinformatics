from pathlib import Path

from mecfs_bio.build_system.wf.base_wf import SimpleWF
from mecfs_bio.build_system.wf.synapse_downloader import SynapseDownloader


class _RecordingDownloader(SynapseDownloader):
    def __init__(self) -> None:
        self.calls: list[tuple[str, Path, str]] = []

    def download(self, synid: str, dest_dir: Path, expected_name: str) -> Path:
        self.calls.append((synid, dest_dir, expected_name))
        out = dest_dir / expected_name
        out.write_text("payload")
        return out


def test_wf_delegates_download_from_synapse(tmp_path: Path):
    fake = _RecordingDownloader()
    wf = SimpleWF(synapse_downloader=fake)

    out = wf.download_from_synapse("syn123", tmp_path, "protein.tar")

    assert out == tmp_path / "protein.tar"
    assert out.read_text() == "payload"
    assert fake.calls == [("syn123", tmp_path, "protein.tar")]
