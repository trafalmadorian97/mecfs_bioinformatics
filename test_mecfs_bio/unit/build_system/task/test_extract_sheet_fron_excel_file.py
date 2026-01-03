from pathlib import Path

import pandas as pd

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.extract_sheet_from_excel_file_task import (
    ExtractSheetFromExelFileTask,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat
from mecfs_bio.build_system.wf.base_wf import SimpleWF


def test_extract_sheet_from_excel_file_task(tmp_path: Path):
    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir()
    source_excel_file_path = tmp_path / "my_excel_file.xlsx"
    sheet_name = "my_sheet"
    dummy_data = pd.DataFrame(
        {
            "a": [1, 2, 3],
            "b": [4, 5, 6],
        }
    )
    dummy_data.to_excel(source_excel_file_path, sheet_name=sheet_name, index=False)
    tsk = ExtractSheetFromExelFileTask(
        meta=SimpleFileMeta("my_extracted_sheet"),
        excel_file_task=FakeTask(meta=SimpleFileMeta("my_fake_task")),
        sheet_name=sheet_name,
        out_format=ParquetOutFormat(),
    )

    def fetch(asset_id: AssetId) -> Asset:
        return FileAsset(source_excel_file_path)

    result = tsk.execute(scratch_dir=scratch_dir, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, FileAsset)
    result_df = pd.read_parquet(result.path)
    pd.testing.assert_frame_equal(dummy_data, result_df)
