import zipfile
from pathlib import Path

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.task.extract_all_zips_in_directory_task import (
    ExtractAllZipsInDir,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import SimpleWF


def test_extract_all_zips(tmp_path: Path):
    name_in_archive = "my_compressed_file"
    name_in_archive_2 = "my_compressed_file_2"
    source_dir = tmp_path / "source_dir"
    source_dir.mkdir(parents=True, exist_ok=True)
    source_zip_path_1 = source_dir / "source_zip_1.zip"
    source_zip_path_2 = source_dir / "source_zip_2.zip"
    scratch_dir = tmp_path / "scratch"
    with zipfile.ZipFile(source_zip_path_1, "w") as source_zip:
        source_zip.writestr(name_in_archive, "12345678910")
    with zipfile.ZipFile(source_zip_path_2, "w") as source_zip:
        source_zip.writestr(name_in_archive_2, "abcdefghij")
    tsk = ExtractAllZipsInDir(
        meta=SimpleDirectoryMeta(AssetId("my_dir")),
        source_directory_task=FakeTask(meta=SimpleDirectoryMeta(AssetId("dummy_dir"))),
    )

    def fetch(asset_id: AssetId) -> Asset:
        return DirectoryAsset(source_dir)

    result = tsk.execute(scratch_dir=scratch_dir, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, DirectoryAsset)
    assert len(list(result.path.iterdir())) == 2
