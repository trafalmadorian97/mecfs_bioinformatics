from pathlib import Path, PurePath

import pandas as pd
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.magma.prepare_gene_sets_for_magma_task import (
    PrepareGeneSetsForMagmaTask,
)
from mecfs_bio.build_system.wf.base_wf import SimpleWF
from mecfs_bio.constants.magma_constants import MAGMA_GENE_COL
from mecfs_bio.constants.msigdb_columns import (
    EXACT_SOURCE,
    NCBI_IDS,
    STANDARD_NAME,
    SYSTEMATIC_NAME,
)
from mecfs_bio.constants.vocabulary_classes.gene_set import MSigDBGeneSetSpec


def _spec(name: str, exact_source: str | None = None) -> MSigDBGeneSetSpec:
    return MSigDBGeneSetSpec(standard_name=name, exact_source=exact_source)


def _write_parquet(tmp_path: Path, rows: list[dict]) -> Path:
    p = tmp_path / "msigdb.parquet"
    pd.DataFrame(rows).to_parquet(p)
    return p


def _make_db_task() -> FakeTask:
    return FakeTask(
        meta=ReferenceFileMeta(
            group="g",
            sub_group="sg",
            sub_folder=PurePath("raw"),
            id=AssetId("db"),
            extension=".parquet",
        )
    )


class _FileFetch(Fetch):
    def __init__(self, path: Path):
        self._path = path

    def __call__(self, asset_id: AssetId) -> Asset:
        return FileAsset(self._path)


_ROWS = [
    {
        STANDARD_NAME: "SET_A",
        SYSTEMATIC_NAME: "M1",
        EXACT_SOURCE: None,
        NCBI_IDS: [1, 2, 3],
    },
    {
        STANDARD_NAME: "SET_B",
        SYSTEMATIC_NAME: "M2",
        EXACT_SOURCE: "SRC:B",
        NCBI_IDS: [2, 3, 4],
    },
    {
        STANDARD_NAME: "SET_C",
        SYSTEMATIC_NAME: "M3",
        EXACT_SOURCE: "SRC:C",
        NCBI_IDS: [5, 6],
    },
]


def test_output_is_tab_separated_with_gene_column(tmp_path: Path):
    db_path = _write_parquet(tmp_path, _ROWS)
    task = PrepareGeneSetsForMagmaTask.create(
        asset_id="out",
        gene_sets=[_spec("SET_A"), _spec("SET_B")],
        parquet_db_task=_make_db_task(),
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch, _FileFetch(db_path), SimpleWF())

    assert isinstance(result, FileAsset)
    df = pd.read_csv(result.path, sep="\t")
    assert MAGMA_GENE_COL in df.columns
    assert "SET_A" in df.columns
    assert "SET_B" in df.columns


def test_gene_column_contains_entrez_ids(tmp_path: Path):
    db_path = _write_parquet(tmp_path, _ROWS)
    task = PrepareGeneSetsForMagmaTask.create(
        asset_id="out",
        gene_sets=[_spec("SET_A"), _spec("SET_B")],
        parquet_db_task=_make_db_task(),
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch, _FileFetch(db_path), SimpleWF())

    assert isinstance(result, FileAsset)
    df = pd.read_csv(result.path, sep="\t")
    # Full universe from _ROWS is {1,2,3,4,5,6}; SET_C genes appear with 0s
    assert set(df[MAGMA_GENE_COL]) == {1, 2, 3, 4, 5, 6}


def test_binary_indicators_are_correct(tmp_path: Path):
    db_path = _write_parquet(tmp_path, _ROWS)
    task = PrepareGeneSetsForMagmaTask.create(
        asset_id="out",
        gene_sets=[_spec("SET_A"), _spec("SET_B")],
        parquet_db_task=_make_db_task(),
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch, _FileFetch(db_path), SimpleWF())

    assert isinstance(result, FileAsset)
    df = pd.read_csv(result.path, sep="\t").set_index(MAGMA_GENE_COL)
    # Gene 1 is only in SET_A
    assert df.loc[1, "SET_A"] == 1
    assert df.loc[1, "SET_B"] == 0
    # Gene 4 is only in SET_B
    assert df.loc[4, "SET_A"] == 0
    assert df.loc[4, "SET_B"] == 1
    # Gene 2 is in both
    assert df.loc[2, "SET_A"] == 1
    assert df.loc[2, "SET_B"] == 1


def test_genes_are_union_across_all_sets(tmp_path: Path):
    db_path = _write_parquet(tmp_path, _ROWS)
    task = PrepareGeneSetsForMagmaTask.create(
        asset_id="out",
        gene_sets=[_spec("SET_A"), _spec("SET_B"), _spec("SET_C", "SRC:C")],
        parquet_db_task=_make_db_task(),
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch, _FileFetch(db_path), SimpleWF())

    assert isinstance(result, FileAsset)
    df = pd.read_csv(result.path, sep="\t")
    assert set(df[MAGMA_GENE_COL]) == {1, 2, 3, 4, 5, 6}


def test_values_are_only_zero_or_one(tmp_path: Path):
    db_path = _write_parquet(tmp_path, _ROWS)
    task = PrepareGeneSetsForMagmaTask.create(
        asset_id="out",
        gene_sets=[_spec("SET_A"), _spec("SET_B")],
        parquet_db_task=_make_db_task(),
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch, _FileFetch(db_path), SimpleWF())

    assert isinstance(result, FileAsset)
    df = pd.read_csv(result.path, sep="\t")
    for col in ["SET_A", "SET_B"]:
        assert set(df[col].unique()).issubset({0, 1})


def test_meta_inherits_reference_file_meta(tmp_path: Path):
    task = PrepareGeneSetsForMagmaTask.create(
        asset_id="my_gene_sets",
        gene_sets=[_spec("SET_A")],
        parquet_db_task=_make_db_task(),
    )
    assert isinstance(task.meta, ReferenceFileMeta)
    assert task.meta.asset_id == AssetId("my_gene_sets")
    assert task.meta.extension == ".txt"
    db_meta = _make_db_task().meta
    assert isinstance(db_meta, ReferenceFileMeta)
    assert task.meta.group == db_meta.group
    assert task.meta.sub_group == db_meta.sub_group
    assert task.meta.sub_folder == db_meta.sub_folder


def test_meta_simple_file_meta_when_dep_is_simple(tmp_path: Path):
    db_task = FakeTask(meta=SimpleFileMeta(AssetId("db")))
    task = PrepareGeneSetsForMagmaTask.create(
        asset_id="my_gene_sets",
        gene_sets=[_spec("SET_A")],
        parquet_db_task=db_task,
    )
    assert isinstance(task.meta, SimpleFileMeta)
    assert task.meta.asset_id == AssetId("my_gene_sets")


def test_deps_is_parquet_db_task(tmp_path: Path):
    db_task = _make_db_task()
    task = PrepareGeneSetsForMagmaTask.create(
        asset_id="out", gene_sets=[_spec("SET_A")], parquet_db_task=db_task
    )
    assert task.deps == [db_task]


def test_unknown_gene_set_raises(tmp_path: Path):
    db_path = _write_parquet(tmp_path, _ROWS)
    task = PrepareGeneSetsForMagmaTask.create(
        asset_id="out",
        gene_sets=[_spec("NONEXISTENT")],
        parquet_db_task=_make_db_task(),
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    with pytest.raises(ValueError, match="NONEXISTENT"):
        task.execute(scratch, _FileFetch(db_path), SimpleWF())


def test_none_ncbi_ids_are_dropped(tmp_path: Path):
    rows_only_null = [
        {
            STANDARD_NAME: "SET_NULL",
            SYSTEMATIC_NAME: "M9",
            EXACT_SOURCE: None,
            NCBI_IDS: [None, 99],
        }
    ]
    db_path = _write_parquet(tmp_path, rows_only_null)
    task = PrepareGeneSetsForMagmaTask.create(
        asset_id="out", gene_sets=[_spec("SET_NULL")], parquet_db_task=_make_db_task()
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch, _FileFetch(db_path), SimpleWF())

    assert isinstance(result, FileAsset)
    df = pd.read_csv(result.path, sep="\t")
    assert set(df[MAGMA_GENE_COL]) == {99}


def test_genes_outside_selected_sets_are_included(tmp_path: Path):
    db_path = _write_parquet(tmp_path, _ROWS)
    task = PrepareGeneSetsForMagmaTask.create(
        asset_id="out", gene_sets=[_spec("SET_A")], parquet_db_task=_make_db_task()
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch, _FileFetch(db_path), SimpleWF())

    assert isinstance(result, FileAsset)
    df = pd.read_csv(result.path, sep="\t").set_index(MAGMA_GENE_COL)
    # Genes 5,6 come from SET_C (not selected) — must be present with indicator 0
    assert 5 in df.index
    assert 6 in df.index
    assert df.loc[5, "SET_A"] == 0
    assert df.loc[6, "SET_A"] == 0
