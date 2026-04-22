"""
Tests implemented by Claude
"""
import sqlite3
from pathlib import Path, PurePath

import pandas as pd
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.sqlite_to_parquet_task import SqliteToParquetTask
from mecfs_bio.build_system.wf.base_wf import SimpleWF


def _make_sqlite_db(tmp_path: Path) -> Path:
    db_path = tmp_path / "test.db"
    conn = sqlite3.connect(db_path)
    conn.executescript("""
        CREATE TABLE gene_set (id INTEGER PRIMARY KEY, name TEXT, collection TEXT);
        CREATE TABLE gene_symbol (id INTEGER PRIMARY KEY, symbol TEXT, ncbi_id TEXT);
        CREATE TABLE gene_set_gene_symbol (gene_set_id INTEGER, gene_symbol_id INTEGER);

        INSERT INTO gene_set VALUES (1, 'SET_A', 'C1');
        INSERT INTO gene_set VALUES (2, 'SET_B', 'C2');
        INSERT INTO gene_symbol VALUES (1, 'BRCA1', '672');
        INSERT INTO gene_symbol VALUES (2, 'TP53', '7157');
        INSERT INTO gene_symbol VALUES (3, 'EGFR', '1956');
        INSERT INTO gene_set_gene_symbol VALUES (1, 1);
        INSERT INTO gene_set_gene_symbol VALUES (1, 2);
        INSERT INTO gene_set_gene_symbol VALUES (2, 3);
    """)
    conn.commit()
    conn.close()
    return db_path


def _make_source_task() -> FakeTask:
    return FakeTask(
        meta=ReferenceFileMeta(
            group="test_group",
            sub_group="test_sub",
            sub_folder=PurePath("raw"),
            id=AssetId("test_sqlite"),
            extension="",
        )
    )


@frozen
class _FileFetch(Fetch):
    path: Path

    def __call__(self, asset_id: AssetId) -> Asset:
        return FileAsset(self.path)


def test_sqlite_to_parquet_simple_query(tmp_path: Path):
    db_path = _make_sqlite_db(tmp_path)
    source_task = _make_source_task()
    scratch = tmp_path / "scratch"
    scratch.mkdir()

    task = SqliteToParquetTask.create(
        source_task=source_task,
        asset_id="out",
        query="SELECT id, name, collection FROM _src.gene_set",
    )
    result = task.execute(scratch_dir=scratch, fetch=_FileFetch(db_path), wf=SimpleWF())

    assert isinstance(result, FileAsset)
    df = pd.read_parquet(result.path)
    assert list(df.columns) == ["id", "name", "collection"]
    assert len(df) == 2
    assert set(df["name"]) == {"SET_A", "SET_B"}


def test_sqlite_to_parquet_with_list_aggregation(tmp_path: Path):
    db_path = _make_sqlite_db(tmp_path)
    source_task = _make_source_task()
    scratch = tmp_path / "scratch"
    scratch.mkdir()

    task = SqliteToParquetTask.create(
        source_task=source_task,
        asset_id="out",
        query="""\
            SELECT
                gs.id,
                gs.name,
                list(sym.symbol ORDER BY gsgs.gene_symbol_id) AS gene_symbols,
                list(TRY_CAST(sym.ncbi_id AS INTEGER) ORDER BY gsgs.gene_symbol_id) AS ncbi_ids
            FROM _src.gene_set gs
            JOIN _src.gene_set_gene_symbol gsgs ON gsgs.gene_set_id = gs.id
            JOIN _src.gene_symbol sym ON sym.id = gsgs.gene_symbol_id
            GROUP BY gs.id, gs.name\
        """,
    )
    result = task.execute(scratch_dir=scratch, fetch=_FileFetch(db_path), wf=SimpleWF())
    assert isinstance(result, FileAsset)

    df = pd.read_parquet(result.path)
    assert "gene_symbols" in df.columns
    assert "ncbi_ids" in df.columns
    row_a = df[df["name"] == "SET_A"].iloc[0]
    assert list(row_a["gene_symbols"]) == ["BRCA1", "TP53"]
    assert list(row_a["ncbi_ids"]) == [672, 7157]
    row_b = df[df["name"] == "SET_B"].iloc[0]
    assert list(row_b["gene_symbols"]) == ["EGFR"]


def test_sqlite_to_parquet_task_meta(tmp_path: Path):
    source_task = _make_source_task()
    task = SqliteToParquetTask.create(
        source_task=source_task,
        asset_id="out_parquet",
        query="SELECT 1 AS x",
    )
    assert isinstance(task.meta, ReferenceFileMeta)
    assert task.meta.extension == ".parquet.zstd"
    assert task.meta.asset_id == AssetId("out_parquet")


def test_sqlite_to_parquet_task_deps(tmp_path: Path):
    source_task = _make_source_task()
    task = SqliteToParquetTask.create(
        source_task=source_task,
        asset_id="out_parquet",
        query="SELECT 1 AS x",
    )
    assert task.deps == [source_task]
