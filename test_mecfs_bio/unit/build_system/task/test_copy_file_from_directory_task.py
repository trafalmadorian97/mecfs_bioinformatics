from pathlib import Path, PurePath

import pandas as pd

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.copy_file_from_directory_task import (
    CopyFileFromDirectoryTask,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import SimpleWF


def test_copy_file_from_directory_task(tmp_path: Path):
    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir()
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    path_inside_directory = PurePath("my_file.csv")
    source_file = source_dir / path_inside_directory
    dummy_data = pd.DataFrame(
        {
            "a": [1, 2, 3],
            "b": [4, 5, 6],
        }
    )
    dummy_data.to_csv(source_file, index=False)
    tsk = CopyFileFromDirectoryTask(
        meta=SimpleFileMeta(
            "myfile",
        ),
        source_directory_task=FakeTask(SimpleDirectoryMeta("my_fake_task")),
        path_inside_directory=PurePath(path_inside_directory),
    )

    def fetch(asset_id: AssetId) -> Asset:
        return DirectoryAsset(source_dir)

    result = tsk.execute(scratch_dir=scratch_dir, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, FileAsset)
