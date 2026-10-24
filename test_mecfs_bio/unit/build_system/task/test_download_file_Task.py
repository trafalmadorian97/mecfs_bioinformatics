import shutil
from pathlib import Path

import pandas as pd
import pytest
from attrs import Factory, frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.download_file_task import (
    DownloadFileTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.build_system.wf.wf_downloader import WFDownloader
from mecfs_bio.util.download.verify import calc_md5_checksum, head_file, verify_hash


@frozen(slots=True)
class FakeWFDownloader(WFDownloader):
    source_file: Path
    requested_connections_log: list[int | None] = Factory(list)

    def download(
        self,
        url: str,
        local_path: Path,
        md5_hash: str | None,
        requested_connections: int | None,
    ) -> None:
        self.requested_connections_log.append(requested_connections)
        shutil.copyfile(self.source_file, local_path)
        verify_hash(
            downloaded_file=local_path,
            expected_hash=md5_hash,
        )


def test_md5_check_passes_fails_appropriately(tmp_path: Path) -> None:
    """
    Check that a mismatched md5 hash raises an error, while a matching hash does not
    """
    external_file_path = tmp_path / "external" / "myfile.txt"
    external_file_path.parent.mkdir(parents=True, exist_ok=True)
    external_file_path.write_text("yoyoyoy")
    hash_of_file = calc_md5_checksum(external_file_path)
    wf = make_wf(downloader=FakeWFDownloader(source_file=external_file_path))
    scratch = tmp_path / "scratch"
    scratch.mkdir(parents=True, exist_ok=True)
    tsk_1 = DownloadFileTask(
        meta=SimpleFileMeta.create("my_file"), url="fake_url", md5_hash=hash_of_file
    )

    tsk_2 = DownloadFileTask(
        meta=SimpleFileMeta.create("my_file"), url="fake_url", md5_hash="wrong_hash"
    )

    def fake_fetch(asset_id: AssetId) -> Asset:
        raise NotImplementedError("")

    tsk_1.execute(scratch_dir=scratch, fetch=fake_fetch, wf=wf)
    with pytest.raises(AssertionError):
        tsk_2.execute(scratch_dir=scratch, fetch=fake_fetch, wf=wf)


def test_requested_connections_reach_downloader(tmp_path: Path) -> None:
    """
    Check that the connections requested by the task are passed through to the downloader
    """
    external_file_path = tmp_path / "external" / "myfile.txt"
    external_file_path.parent.mkdir(parents=True, exist_ok=True)
    external_file_path.write_text("yoyoyoy")
    downloader = FakeWFDownloader(source_file=external_file_path)
    scratch = tmp_path / "scratch"
    scratch.mkdir(parents=True, exist_ok=True)
    num_connections = 8
    tsk = DownloadFileTask(
        meta=SimpleFileMeta.create("my_file"),
        url="fake_url",
        md5_hash=None,
        connections=num_connections,
    )

    def fake_fetch(asset_id: AssetId) -> Asset:
        raise NotImplementedError("")

    tsk.execute(
        scratch_dir=scratch, fetch=fake_fetch, wf=make_wf(downloader=downloader)
    )
    assert downloader.requested_connections_log == [num_connections]


def test_head(tmp_path: Path) -> None:
    """
    Verify that we do not get an error when attempting to print the head of a non text file.
    """
    dummy_parquet = tmp_path / "dummy_parquet.parquet"
    pd.DataFrame({"a": [1, 2, 3]}).to_parquet(dummy_parquet)
    head_file(dummy_parquet)
