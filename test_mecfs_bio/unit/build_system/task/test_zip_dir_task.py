import shutil
from pathlib import Path

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.zip_dir_task import ZipDirTask
from mecfs_bio.build_system.wf.base_wf import make_wf


def test_zip_dir_task(tmp_path: Path):
    to_zip = tmp_path / "to_Zip"
    to_zip.mkdir()
    file_1 = to_zip / "file_1.txt"
    file_2 = to_zip / "file_2.txt"
    file_1.write_text("hello")
    file_2.write_text("goodbye")
    unpacked_result = tmp_path / "unpacked_result"
    tsk = ZipDirTask(
        source_dir_task=FakeTask(SimpleDirectoryMeta(AssetId("source_dir")), deps=[]),
        meta=SimpleFileMeta(AssetId("zipped")),
    )
    scratch_dir = to_zip / "scratch_dir"

    def fetch(asset_id: AssetId) -> Asset:
        return DirectoryAsset(to_zip)

    result = tsk.execute(scratch_dir=scratch_dir, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    shutil.unpack_archive(result.path, unpacked_result)
    assert (unpacked_result / "file_1.txt").read_text() == "hello"
    assert (unpacked_result / "file_2.txt").read_text() == "goodbye"
