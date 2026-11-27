import gzip
import math

import numpy as np
import polars as pl
import pyarrow.parquet as pq
import pytest

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.decode_database.build_decode_somascan_protein_task import (
    DECODE_POSMAP_FILENAME,
    DECODE_SEEK_INDEX_FILENAME,
    DECODE_SIGNIFICANT_ROWS_FILENAME,
    DECODE_SLIM_FILENAME,
    DECODE_SLIM_MEMBER,
    SIGNIFICANT_ROWS_SCANNED_KEY,
    SLIM_INDEX_ROWS_MISSING_KEY,
)
from mecfs_bio.build_system.task.decode_database.decode_sumstats_format import (
    chrom_code,
)
from mecfs_bio.build_system.task.decode_database.one_pass import (
    POSMAP_CHROM_COL,
    POSMAP_LINE_OFFSET_COL,
    POSMAP_POS_COL,
    POSMAP_ROW_NUMBER_COL,
    SIGNIFICANT_ROW_NUMBER_COL,
    SIGNIFICANT_ROWS_COLUMNS,
)
from mecfs_bio.build_system.task.directory_member_dataframe import (
    scan_directory_member,
)
from mecfs_bio.constants.decode_somascan_constants import (
    DECODE_INDEX_NAME_COL,
    DECODE_MLOG10P_COL,
    DECODE_SIGNIFICANT_ROWS_MIN_MLOG10P,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)
from mecfs_bio.util.file_io.gzip_seek_index import seek_points
from test_mecfs_bio.unit.build_system.task.decode_database.decode_fixtures import (
    DecodeDatabaseFixture,
    decode_database_fixture,
    decode_text,
)

_PROTEIN_KEY = "final_somascan_smp/protein_a.txt.gz"


@pytest.fixture(scope="module")
def database(tmp_path_factory: pytest.TempPathFactory) -> DecodeDatabaseFixture:
    return decode_database_fixture(tmp_path_factory.mktemp("decode"), seed=11)


@pytest.fixture(scope="module")
def built(database: DecodeDatabaseFixture) -> DirectoryAsset:
    return database.build_protein(database.protein_task(database.rows, _PROTEIN_KEY))


def _index(database: DecodeDatabaseFixture) -> pl.DataFrame:
    asset = database.assets[database.index_task.asset_id]
    assert isinstance(asset, FileAsset)
    return pl.read_parquet(asset.path)


def _expected_slim(index: pl.DataFrame, protein_rows: pl.DataFrame) -> pl.DataFrame:
    """Beta/SE/N of the protein's row with each index NAME, null where it has none."""
    return (
        index.select(DECODE_INDEX_NAME_COL)
        .join(
            protein_rows,
            left_on=DECODE_INDEX_NAME_COL,
            right_on="Name",
            how="left",
            maintain_order="left",
        )
        .select(
            pl.col("Beta").cast(pl.Float32).alias(GWASLAB_BETA_COL),
            pl.col("SE").cast(pl.Float32).alias(GWASLAB_SE_COL),
            pl.col("N").cast(pl.Float32).alias(GWASLAB_SAMPLE_SIZE_COLUMN),
        )
    )


def _slim_missing(built: DirectoryAsset) -> int:
    metadata = pq.read_schema(built.path / DECODE_SLIM_FILENAME).metadata
    assert metadata is not None
    return int(metadata[SLIM_INDEX_ROWS_MISSING_KEY.encode()])


def test_protein_outputs_match_brute_force(
    database: DecodeDatabaseFixture, built: DirectoryAsset
) -> None:
    rows = database.rows
    assert pl.read_parquet(built.path / DECODE_SLIM_FILENAME).equals(
        _expected_slim(_index(database), rows)
    )
    assert _slim_missing(built) == 0

    text = decode_text(rows)
    posmap = pl.read_parquet(built.path / DECODE_POSMAP_FILENAME)
    assert posmap.height > 50
    for row in posmap.iter_rows(named=True):
        offset = row[POSMAP_LINE_OFFSET_COL]
        assert text[offset - 1 : offset] == b"\n"
        fields = text[offset : offset + 64].split(b"\t")
        assert chrom_code(fields[0].decode()) == row[POSMAP_CHROM_COL]
        assert int(fields[1]) == row[POSMAP_POS_COL]
        assert (
            rows.row(row[POSMAP_ROW_NUMBER_COL], named=True)["Pos"]
            == row[POSMAP_POS_COL]
        )

    source = database.objects[
        next(r for r in database.objects if r.key == _PROTEIN_KEY)
    ]
    with open(source, "rb") as f:
        points = seek_points((built.path / DECODE_SEEK_INDEX_FILENAME).read_bytes(), f)
    assert len(points) > 10
    assert gzip.decompress(source.read_bytes()) == text


def test_significant_rows_match_brute_force(
    database: DecodeDatabaseFixture, built: DirectoryAsset
) -> None:
    rows = database.rows.with_row_index(SIGNIFICANT_ROW_NUMBER_COL)
    expected = rows.filter(
        pl.col(DECODE_MLOG10P_COL) >= DECODE_SIGNIFICANT_ROWS_MIN_MLOG10P
    )
    path = built.path / DECODE_SIGNIFICANT_ROWS_FILENAME
    got = pl.read_parquet(path)
    assert got.columns == SIGNIFICANT_ROWS_COLUMNS
    assert got.equals(
        expected.select(SIGNIFICANT_ROWS_COLUMNS).with_columns(
            pl.col(SIGNIFICANT_ROW_NUMBER_COL).cast(pl.UInt64)
        )
    )
    assert set(got["Chrom"]) >= {"chrX", "chr2"}
    probes = got.filter(pl.col("Chrom") == "chrX")[DECODE_MLOG10P_COL].to_list()
    assert any(math.isclose(p, 7.31) for p in probes)
    assert not any(math.isclose(p, 7.29) for p in probes)
    metadata = pq.read_schema(path).metadata
    assert int(metadata[SIGNIFICANT_ROWS_SCANNED_KEY.encode()]) == database.rows.height


def _index_row_numbers(database: DecodeDatabaseFixture) -> np.ndarray:
    names = _index(database)[DECODE_INDEX_NAME_COL]
    return np.flatnonzero(database.rows["Name"].is_in(names.implode()).to_numpy())


def test_protein_with_reordered_and_deleted_rows(
    database: DecodeDatabaseFixture,
) -> None:
    """
    Rows outside the index deleted and two index rows swapped: each index row still gets its
    own row's values.
    """
    in_index = _index_row_numbers(database)
    outside = np.setdiff1d(np.arange(database.rows.height), in_index)
    order = np.setdiff1d(np.arange(database.rows.height), outside[[5, 40, 41, 300]])
    i, j = (int(np.flatnonzero(order == r)[0]) for r in in_index[[10, 11]])
    order[i], order[j] = order[j], order[i]
    protein_rows = database.rows[order]
    built = database.build_protein(
        database.protein_task(protein_rows, "final_somascan_smp/reordered.txt.gz")
    )
    assert pl.read_parquet(built.path / DECODE_SLIM_FILENAME).equals(
        _expected_slim(_index(database), protein_rows)
    )
    assert _slim_missing(built) == 0


def test_missing_index_rows_are_null(database: DecodeDatabaseFixture) -> None:
    in_index = _index_row_numbers(database)
    dropped = in_index[[3, 200]]
    protein_rows = database.rows[np.setdiff1d(np.arange(database.rows.height), dropped)]
    built = database.build_protein(
        database.protein_task(protein_rows, "final_somascan_smp/missing.txt.gz")
    )
    slim = pl.read_parquet(built.path / DECODE_SLIM_FILENAME)
    assert slim.equals(_expected_slim(_index(database), protein_rows))
    assert slim[GWASLAB_BETA_COL].null_count() == 2
    assert _slim_missing(built) == 2


def test_protein_fails_on_duplicated_index_name(
    database: DecodeDatabaseFixture,
) -> None:
    in_index = _index_row_numbers(database)
    order = np.sort(np.concatenate([np.arange(database.rows.height), in_index[[7]]]))
    task = database.protein_task(
        database.rows[order], "final_somascan_smp/duplicated.txt.gz"
    )
    with pytest.raises(pl.exceptions.ComputeError):
        database.build_protein(task)


def test_protein_fails_when_name_alleles_disagree(
    database: DecodeDatabaseFixture,
) -> None:
    """A row keeping an index Name but with its alleles swapped must not be copied."""
    row = int(_index_row_numbers(database)[9])
    swapped = database.rows.with_columns(
        pl.when(pl.int_range(pl.len()) == row)
        .then(pl.col(col))
        .otherwise(pl.col(other))
        .alias(other)
        for col, other in (
            ("effectAllele", "otherAllele"),
            ("otherAllele", "effectAllele"),
        )
    )
    task = database.protein_task(swapped, "final_somascan_smp/swapped_alleles.txt.gz")
    with pytest.raises(AssertionError):
        database.build_protein(task)


def test_protein_fails_on_md5_mismatch(database: DecodeDatabaseFixture) -> None:
    task = database.protein_task(
        database.rows, "final_somascan_smp/stale.txt.gz", md5="0" * 32
    )
    with pytest.raises(AssertionError):
        database.build_protein(task)


def test_slim_reads_through_directory_member(built: DirectoryAsset) -> None:
    slim = scan_directory_member(built, DECODE_SLIM_MEMBER).to_native().collect()
    assert slim.columns == [
        GWASLAB_BETA_COL,
        GWASLAB_SE_COL,
        GWASLAB_SAMPLE_SIZE_COLUMN,
    ]
