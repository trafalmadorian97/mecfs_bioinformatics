from pathlib import Path

import pandas as pd

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.wf.base_wf import SimpleWF


def test_join_dataframes_task(tmp_path: Path):
    df_1_loc = tmp_path / "df_1.csv"
    df_2_loc = tmp_path / "df_2.csv"
    df_1 = pd.DataFrame({"col1": [1, 2, 3], "col2": [4, 5, 6]})
    df_1.to_csv(df_1_loc, index=False)
    df_2 = pd.DataFrame({"col1": [1, 2, 3], "col3": ["a", "b", "c"]})
    df_2.to_csv(df_2_loc, index=False)
    scratch_loc = tmp_path / "scratch"
    scratch_loc.mkdir(exist_ok=True, parents=True)
    task = JoinDataFramesTask(
        df_1_task=FakeTask(
            SimpleFileMeta(
                AssetId("df_1_task"),
                read_spec=DataFrameReadSpec(DataFrameTextFormat(",")),
            )
        ),
        df_2_task=FakeTask(
            SimpleFileMeta(
                AssetId("df_2_task"),
                read_spec=DataFrameReadSpec(DataFrameTextFormat(",")),
            )
        ),
        how="inner",
        left_on=["col1"],
        right_on=["col1"],
        meta=SimpleFileMeta(AssetId("joined")),
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "df_1_task":
            return FileAsset(df_1_loc)
        if asset_id == "df_2_task":
            return FileAsset(df_2_loc)
        raise ValueError("unknown asset id")

    result = task.execute(scratch_dir=scratch_loc, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, FileAsset)
    result_df = pd.read_csv(result.path)
    assert len(result_df) == 3
