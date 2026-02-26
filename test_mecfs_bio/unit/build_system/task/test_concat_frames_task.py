from pathlib import Path
import pandas as pd

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec, DataFrameTextFormat
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.concat_frames_task import ConcatFramesTask
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat
from mecfs_bio.build_system.wf.base_wf import SimpleWF


def test_concat_frames_in_dir_task(tmp_path: Path):
    subdir = tmp_path / "subdir"
    subdir.mkdir()
    df_1_loc = subdir / "df_1.csv"
    df_2_loc = subdir / "df_2.csv"
    df_1 = pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "d", "z"]})
    df_1.to_csv(df_1_loc, index=False)
    df_2 = pd.DataFrame({"col1": [1, 2, 3], "col2": ["a", "b", "c"]})
    df_2.to_csv(df_2_loc, index=False)
    scratch_loc = tmp_path / "scratch"

    scratch_loc.mkdir(exist_ok=True, parents=True)
    task = ConcatFramesTask(
        meta=SimpleFileMeta("concat_result"),
        frames_tasks=[
            FakeTask(
                SimpleFileMeta("frame1",
                               read_spec=DataFrameReadSpec(DataFrameTextFormat(",")))
            ),

            FakeTask(
                SimpleFileMeta("frame2",
                               read_spec=DataFrameReadSpec(DataFrameTextFormat(",")))
            ),
        ],
        out_format=ParquetOutFormat()
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "frame1":
            return FileAsset(df_1_loc)
        if asset_id == "frame2":
            return FileAsset(df_2_loc)
        raise ValueError("unknown asset id")

    result = task.execute(scratch_dir=scratch_loc, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, FileAsset)
    result_df = pd.read_parquet(result.path)
    assert len(result_df) == 6