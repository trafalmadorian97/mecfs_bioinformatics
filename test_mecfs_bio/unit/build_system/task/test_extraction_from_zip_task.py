import os
import zipfile
from pathlib import Path

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.extraction_one_file_from_zip_task import (
    ExtractFromZipTask,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.make_executable_wrapper_task import (
    MakeExecutableWrapperTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf


def test_extract_from_zip_and_make_executable(tmp_path: Path):
    name_in_archive = "my_compressed_file"
    source_zip_path = tmp_path / "source_zip.zip"
    scratch_dir = tmp_path / "scratch"
    with zipfile.ZipFile(source_zip_path, "w") as source_zip:
        source_zip.writestr(name_in_archive, "12345678910")
    tsk = MakeExecutableWrapperTask(
        ExtractFromZipTask(
            meta=SimpleFileMeta(AssetId("my_file")),
            source_file_task=FakeTask(meta=SimpleFileMeta(AssetId("dummy_file"))),
            file_to_extract=name_in_archive,
        )
    )

    def fetch(asset_id: AssetId) -> Asset:
        return FileAsset(source_zip_path)

    result = tsk.execute(scratch_dir=scratch_dir, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    assert result.path == scratch_dir / name_in_archive
    assert os.access(result.path, os.X_OK)
