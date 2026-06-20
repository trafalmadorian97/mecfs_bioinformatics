from pathlib import Path

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.copy_files_into_directory_task import (
    CopyFilesIntoDirectoryTask,
    CopySource,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import SimpleWF


def test_copy_files(tmp_path: Path):
    file_1_loc = tmp_path / "file_1.txt"
    file_1_loc.write_text("hello world")
    file_1_task = FakeTask(meta=SimpleFileMeta(AssetId("file1")))
    file_1_in_dir_name = "name1"
    file_2_loc = tmp_path / "file_2.txt"
    file_2_loc.write_text("hello world_2")
    file_2_task = FakeTask(meta=SimpleFileMeta(AssetId("file2")))
    file_2_in_dir_name = "name2"
    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir()
    tsk = CopyFilesIntoDirectoryTask(
        sources=[
            CopySource(file_1_task, name_in_dir=file_1_in_dir_name),
            CopySource(file_2_task, name_in_dir=file_2_in_dir_name),
        ],
        meta=SimpleFileMeta(AssetId("Dir")),
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "file1":
            return FileAsset(file_1_loc)
        if asset_id == "file2":
            return FileAsset(file_2_loc)
        raise ValueError("Invalid asset_id")

    result = tsk.execute(scratch_dir=scratch_dir, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, DirectoryAsset)
    assert (result.path / file_1_in_dir_name).exists()
    assert (result.path / file_2_in_dir_name).exists()
