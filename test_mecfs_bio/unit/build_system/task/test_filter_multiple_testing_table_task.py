from pathlib import Path

import pandas as pd
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.fdr_multiple_testing_table_task import (
    FilterMultipleTestingTableTask,
    Method,
)
from mecfs_bio.build_system.wf.base_wf import SimpleWF

_dummy_df = pd.DataFrame({"p": [0.02, 0.02, 1, 1, 1]})


@pytest.mark.parametrize(
    ["procedure", "expected_passed"], [["fdr_bh", 2], ["bonferroni", 0]]
)
def test_multiple_testing(tmp_path: Path, procedure: Method, expected_passed: int):
    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir(exist_ok=True, parents=True)
    input_path = tmp_path / "input.csv"
    _dummy_df.to_csv(input_path, index=False)
    tsk = FilterMultipleTestingTableTask(
        meta=SimpleFileMeta(
            AssetId("my_filtered_df"),
            read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator=",")),
        ),
        table_source_task=FakeTask(
            SimpleFileMeta(
                AssetId("dummy"),
                read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator=",")),
            ),
        ),
        p_value_column="p",
        alpha=0.05,
        method=procedure,
    )

    def fetch(asset_id: AssetId) -> Asset:
        return FileAsset(input_path)

    result = tsk.execute(scratch_dir=scratch_dir, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, FileAsset)
    result_df = pd.read_csv(result.path)
    assert len(result_df) == expected_passed
