import json
from pathlib import Path, PurePath

import pandas as pd

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.json_to_parquet_task import JsonToParquetTask
from mecfs_bio.build_system.wf.base_wf import SimpleWF


def _make_json_task(tmp_path: Path, data: list[dict]) -> tuple[FakeTask, Path]:
    json_path = tmp_path / "input.json"
    json_path.write_text(json.dumps(data))
    task = FakeTask(
        meta=ReferenceFileMeta(
            group="test_group",
            sub_group="test_sub_group",
            sub_folder=PurePath("raw"),
            id=AssetId("test_json"),
            extension=".json",
        )
    )
    return task, json_path


class _FileFetch(Fetch):
    def __init__(self, path: Path):
        self._path = path

    def __call__(self, asset_id: AssetId) -> Asset:
        return FileAsset(self._path)


def test_json_to_parquet_task_executes(tmp_path: Path):
    records = [{"a": 1, "b": "x"}, {"a": 2, "b": "y"}]
    json_task, json_path = _make_json_task(tmp_path, records)
    scratch = tmp_path / "scratch"
    scratch.mkdir()

    task = JsonToParquetTask.create(json_task=json_task, asset_id="test_parquet")
    result = task.execute(
        scratch_dir=scratch, fetch=_FileFetch(json_path), wf=SimpleWF()
    )

    assert isinstance(result, FileAsset)
    df = pd.read_parquet(result.path)
    assert list(df.columns) == ["a", "b"]
    assert len(df) == 2


def test_json_to_parquet_task_meta(tmp_path: Path):
    records = [{"x": 10}]
    json_task, _ = _make_json_task(tmp_path, records)

    task = JsonToParquetTask.create(json_task=json_task, asset_id="parquet_out")
    assert isinstance(task.meta, ReferenceFileMeta)
    assert task.meta.extension == ".parquet.zstd"
    assert task.meta.asset_id == AssetId("parquet_out")


def test_json_to_parquet_task_deps(tmp_path: Path):
    records = [{"x": 10}]
    json_task, _ = _make_json_task(tmp_path, records)
    task = JsonToParquetTask.create(json_task=json_task, asset_id="parquet_out")
    assert task.deps == [json_task]
