from pathlib import Path

import pytest

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.decode_database.decode_somascan_protein import (
    DecodeSomascanProtein,
)
from mecfs_bio.build_system.task.decode_database.download_decode_somascan_protein_task import (
    DownloadDecodeSomascanProteinTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from test_mecfs_bio.unit.build_system.task.decode_database.decode_fixtures import (
    FakeParallelRangeDownloader,
    FakeS3Access,
    synthetic_protein,
    synthetic_rows,
    write_synthetic_decode_file,
)

_KEY = "final_somascan_smp/protein_d.txt.gz"


def _run(tmp_path: Path, protein: DecodeSomascanProtein, access: FakeS3Access) -> Path:
    task = DownloadDecodeSomascanProteinTask.create("protein_d_raw", protein)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    asset = task.execute(
        scratch,
        fetch=lambda asset_id: pytest.fail(f"unexpected fetch {asset_id}"),
        wf=make_wf(
            s3_access=access, parallel_range_downloader=FakeParallelRangeDownloader()
        ),
    )
    assert isinstance(asset, FileAsset)
    return asset.path


def test_download_keeps_the_verified_file(tmp_path: Path) -> None:
    source = tmp_path / "source.txt.gz"
    write_synthetic_decode_file(source, synthetic_rows(seed=3))
    protein = synthetic_protein(source, _KEY)
    out = _run(tmp_path, protein, FakeS3Access(objects={protein.ref: source}))
    assert out.read_bytes() == source.read_bytes()


def test_download_fails_when_bucket_md5_changed(tmp_path: Path) -> None:
    source = tmp_path / "source.txt.gz"
    write_synthetic_decode_file(source, synthetic_rows(seed=3))
    protein = synthetic_protein(source, _KEY)
    access = FakeS3Access(
        objects={protein.ref: source}, reported_md5={protein.ref: "0" * 32}
    )
    with pytest.raises(AssertionError):
        _run(tmp_path, protein, access)
