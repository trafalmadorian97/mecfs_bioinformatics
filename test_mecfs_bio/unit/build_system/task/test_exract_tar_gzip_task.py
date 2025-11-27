import tarfile
import tempfile
from pathlib import Path
from typing import Literal

import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.extract_tar_gzip_task import (
    ExtractTarGzipTask,
    ReadMode,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import SimpleWF

WriteMode = Literal["w:gz", "w"]


def _create_dummy_tar_file(
    contents: str, name_in_tar: str, pth: Path, mode: WriteMode = "w:gz"
) -> None:
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir) / name_in_tar
        tmp_path.write_text(contents)
        with tarfile.open(pth, mode=mode) as tar_object:
            tar_object.add(tmp_path, arcname=name_in_tar)


@pytest.mark.parametrize(
    ["read_mode", "write_mode"],
    [
        ["r:gz", "w:gz"],
        ["r", "w"],
    ],
)
def test_extract_tar_gzip_task(
    tmp_path: Path, read_mode: ReadMode, write_mode: WriteMode
):
    tar_loc = tmp_path / "my_tar.tar.gzip"
    file_contents = "abc123"
    name_in_tar = "my_tarred_file"
    _create_dummy_tar_file(
        file_contents, pth=tar_loc, name_in_tar=name_in_tar, mode=write_mode
    )
    scratch_dir = tmp_path / "scratch"
    tsk = ExtractTarGzipTask(
        meta=SimpleFileMeta(AssetId("my_file")),
        source_file_task=FakeTask(meta=SimpleFileMeta(AssetId("dummy_file"))),
        subdir_name=None,
        read_mode=read_mode,
    )

    def fetch(asset_id: AssetId) -> Asset:
        return FileAsset(tar_loc)

    result = tsk.execute(scratch_dir=scratch_dir, fetch=fetch, wf=SimpleWF())
    assert result.path == scratch_dir
    dir_contents = list(scratch_dir.glob("*"))
    print(dir_contents)
    assert (scratch_dir / name_in_tar).exists()
    assert (scratch_dir / name_in_tar).read_text() == file_contents
