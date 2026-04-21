from pathlib import Path, PurePath

import duckdb
import pandas as pd

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.unpack_map_parquet_task import UnpackMapParquetTask
from mecfs_bio.build_system.wf.base_wf import SimpleWF


def _make_map_parquet(tmp_path: Path) -> Path:
    """Write a parquet file with a MAP(VARCHAR, STRUCT) column matching the MSigDB shape."""
    parquet_path = tmp_path / "map_input.parquet"
    duckdb.sql(f"""
        COPY (
            SELECT MAP {{
                'gene_set_A': {{'collection': 'C1', 'geneSymbols': ['BRCA1', 'TP53']}},
                'gene_set_B': {{'collection': 'C2', 'geneSymbols': ['EGFR']}}
            }} AS data
        ) TO '{parquet_path}' (FORMAT 'PARQUET')
    """)
    return parquet_path


def _make_source_task() -> FakeTask:
    return FakeTask(
        meta=ReferenceFileMeta(
            group="test_group",
            sub_group="test_sub",
            sub_folder=PurePath("raw"),
            id=AssetId("map_parquet"),
            extension=".parquet.zstd",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
    )


class _FileFetch(Fetch):
    def __init__(self, path: Path):
        self._path = path

    def __call__(self, asset_id: AssetId) -> Asset:
        return FileAsset(self._path)


def test_unpack_map_parquet_task_executes(tmp_path: Path):
    parquet_path = _make_map_parquet(tmp_path)
    source_task = _make_source_task()
    scratch = tmp_path / "scratch"
    scratch.mkdir()

    task = UnpackMapParquetTask.create(
        source_task=source_task,
        asset_id="unpacked",
        map_column="data",
        name_column="name",
    )
    result = task.execute(
        scratch_dir=scratch, fetch=_FileFetch(parquet_path), wf=SimpleWF()
    )

    assert isinstance(result, FileAsset)
    df = pd.read_parquet(result.path)
    assert "name" in df.columns
    assert "collection" in df.columns
    assert "geneSymbols" in df.columns
    assert len(df) == 2
    assert set(df["name"]) == {"gene_set_A", "gene_set_B"}


def test_unpack_map_parquet_task_gene_symbols_are_lists(tmp_path: Path):
    parquet_path = _make_map_parquet(tmp_path)
    source_task = _make_source_task()
    scratch = tmp_path / "scratch"
    scratch.mkdir()

    task = UnpackMapParquetTask.create(
        source_task=source_task,
        asset_id="unpacked",
        map_column="data",
    )
    result = task.execute(
        scratch_dir=scratch, fetch=_FileFetch(parquet_path), wf=SimpleWF()
    )
    assert isinstance(result, FileAsset)
    df = pd.read_parquet(result.path)
    row_a = df[df["name"] == "gene_set_A"].iloc[0]
    assert list(row_a["geneSymbols"]) == ["BRCA1", "TP53"]


def test_unpack_map_parquet_task_meta(tmp_path: Path):
    source_task = _make_source_task()
    task = UnpackMapParquetTask.create(source_task=source_task, asset_id="unpacked_out")
    assert isinstance(task.meta, ReferenceFileMeta)
    assert task.meta.extension == ".parquet.zstd"
    assert task.meta.asset_id == AssetId("unpacked_out")


def test_unpack_map_parquet_task_deps(tmp_path: Path):
    source_task = _make_source_task()
    task = UnpackMapParquetTask.create(source_task=source_task, asset_id="unpacked_out")
    assert task.deps == [source_task]
