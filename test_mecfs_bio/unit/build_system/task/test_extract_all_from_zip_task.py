import zipfile
from pathlib import Path

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.extract_all_from_zip_task import ExtractAllFromZipTask
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import SimpleWF


def test_extract_all_zip(tmp_path: Path):
    name_in_archive = "my_compressed_file"
    name_in_archive_2 = "my_compressed_file_2"
    source_zip_path = tmp_path / "source_zip.zip"
    scratch_dir = tmp_path / "scratch"
    with zipfile.ZipFile(source_zip_path, "w") as source_zip:
        source_zip.writestr(name_in_archive, "12345678910")
        source_zip.writestr(name_in_archive_2, "abcdefghij")
    tsk = ExtractAllFromZipTask(
        meta=SimpleDirectoryMeta(AssetId("my_dir")),
        source_file_task=FakeTask(meta=SimpleFileMeta(AssetId("dummy_file"))),
    )

    def fetch(asset_id: AssetId) -> Asset:
        return FileAsset(source_zip_path)

    result = tsk.execute(scratch_dir=scratch_dir, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, DirectoryAsset)
    assert len(list(result.path.iterdir())) == 2


def test_extract_all_zip_sub_folder(tmp_path: Path):
    inner_folder = "inner"
    source_zip_path = tmp_path / "source_zip.zip"
    scratch_dir = tmp_path / "scratch"
    with zipfile.ZipFile(source_zip_path, "w") as source_zip:
        source_zip.writestr(f"{inner_folder}/file_a", "12345678910")
        source_zip.writestr(f"{inner_folder}/file_b", "abcdefghij")
    tsk = ExtractAllFromZipTask(
        meta=SimpleDirectoryMeta(AssetId("my_dir")),
        source_file_task=FakeTask(meta=SimpleFileMeta(AssetId("dummy_file"))),
        sub_folder_name_inside_zip=inner_folder,
    )

    def fetch(asset_id: AssetId) -> Asset:
        return FileAsset(source_zip_path)

    result = tsk.execute(scratch_dir=scratch_dir, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, DirectoryAsset)
    # The returned directory is the inner folder, holding its two files directly.
    assert result.path.name == inner_folder
    assert {p.name for p in result.path.iterdir()} == {"file_a", "file_b"}
