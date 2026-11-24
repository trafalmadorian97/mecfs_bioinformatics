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
    DECODE_INDEX_SOURCE_ROW_COL,
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


def test_protein_outputs_match_brute_force(
    database: DecodeDatabaseFixture, built: DirectoryAsset
) -> None:
    rows = database.rows
    source_rows = _index(database)[DECODE_INDEX_SOURCE_ROW_COL]
    expected = rows[source_rows].select(
        pl.col("Beta").cast(pl.Float32).alias(GWASLAB_BETA_COL),
        pl.col("SE").cast(pl.Float32).alias(GWASLAB_SE_COL),
        pl.col("N").cast(pl.Float32).alias(GWASLAB_SAMPLE_SIZE_COLUMN),
    )
    assert pl.read_parquet(built.path / DECODE_SLIM_FILENAME).equals(expected)

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


def _swap_two_index_rows(database: DecodeDatabaseFixture) -> pl.DataFrame:
    source_rows = _index(database)[DECODE_INDEX_SOURCE_ROW_COL].to_numpy()
    i, j = int(source_rows[10]), int(source_rows[11])
    order = np.arange(database.rows.height)
    order[i], order[j] = j, i
    return database.rows[order]


def test_protein_fails_on_swapped_rows(database: DecodeDatabaseFixture) -> None:
    task = database.protein_task(
        _swap_two_index_rows(database), "final_somascan_smp/swapped.txt.gz"
    )
    with pytest.raises(AssertionError):
        database.build_protein(task)


def test_protein_fails_on_wrong_total_rows(database: DecodeDatabaseFixture) -> None:
    task = database.protein_task(
        database.rows.head(-1), "final_somascan_smp/short.txt.gz"
    )
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
