from pathlib import Path, PurePath

import pandas as pd
import pytest

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.upset_plot_task import (
    DirSetSource,
    FileSetSource,
    UpSetPlotTask,
)
from mecfs_bio.build_system.wf.base_wf import SimpleWF


@pytest.mark.parametrize(
    ["df_1", "df_2"],
    [
        (
            pd.DataFrame({"col1": ["A", "B", "C"]}),
            pd.DataFrame({"col1": ["B", "C", "D"]}),
        ),
        (pd.DataFrame({"col1": []}), pd.DataFrame({"col1": []})),
    ],
)
def test_upset_plot_task(
    tmp_path: Path,
    df_1: pd.DataFrame,
    df_2: pd.DataFrame,
):
    scratch = tmp_path / "scratch"
    scratch.mkdir(parents=True, exist_ok=True)
    path_in_dir = PurePath("file_in_dir.parquet")
    file1_path = tmp_path / "file1.parquet"
    df_1.to_parquet(file1_path)
    dir1_path = tmp_path / "dir"
    dir1_path.mkdir()
    df_2.to_parquet(dir1_path / path_in_dir)
    plot_task = UpSetPlotTask(
        meta=SimpleFileMeta(id="my_aset_id"),
        set_sources=[
            FileSetSource(
                "file_source",
                task=FakeTask(
                    meta=SimpleFileMeta(
                        "file1", read_spec=DataFrameReadSpec(DataFrameParquetFormat())
                    )
                ),
                col_name="col1",
            ),
            DirSetSource(
                name="dir_source",
                task=FakeTask(SimpleDirectoryMeta("dir1")),
                file_in_dir=path_in_dir,
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
                col_name="col1",
            ),
        ],
    )

    def fetch(asset_id: AssetId):
        if asset_id == "file1":
            return FileAsset(file1_path)
        if asset_id == "dir1":
            return DirectoryAsset(dir1_path)
        raise ValueError(f"Unknown asset id {asset_id}")

    result = plot_task.execute(scratch_dir=scratch, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, FileAsset)
