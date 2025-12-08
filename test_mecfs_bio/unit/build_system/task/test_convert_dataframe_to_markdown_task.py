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
from mecfs_bio.build_system.task.convert_dataframe_to_markdown_task import (
    ConvertDataFrameToMarkdownTask,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import SimpleWF


def test_convert_dataframe_to_markdown_task(tmp_path: Path):
    df_1_loc = tmp_path / "df_1.csv"
    df_1 = pd.DataFrame({"col1": [1, 2, 3], "col2": [4, 5, 6]})
    df_1.to_csv(df_1_loc, index=False)
    scratch_loc = tmp_path / "scratch"
    scratch_loc.mkdir(exist_ok=True, parents=True)
    task = ConvertDataFrameToMarkdownTask(
        meta=SimpleFileMeta(AssetId("my_task")),
        df_task=FakeTask(
            SimpleFileMeta(
                AssetId("Fake"),
                read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
            ),
        ),
    )

    def fetch(asset_id: AssetId) -> Asset:
        return FileAsset(df_1_loc)

    result = task.execute(scratch_dir=scratch_loc, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, FileAsset)
    assert result.path.is_file()
