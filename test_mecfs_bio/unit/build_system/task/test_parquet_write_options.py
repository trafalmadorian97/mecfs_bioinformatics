from pathlib import Path

import polars as pl
import pyarrow.parquet as pq

from mecfs_bio.build_system.task.pipe_dataframe_task import (
    ParquetWriteOptions,
    write_parquet_with_options,
)


def _frame() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "label": ["a", "b", "a", "b"],
            "count": [1, 2, 3, 4],
            "value": [0.125, 0.25, 0.5, 1.0],
        }
    )


def _table():
    return _frame().to_arrow()


def _float_encodings(path: Path) -> set[str]:
    row_group = pq.ParquetFile(path).metadata.row_group(0)
    columns = [row_group.column(i) for i in range(row_group.num_columns)]
    value_column = next(c for c in columns if c.path_in_schema == "value")
    return set(value_column.encodings)


def test_byte_stream_split_is_applied_to_float_columns(tmp_path: Path):
    # Dictionary encoding silently wins over BYTE_STREAM_SPLIT if it is left
    # enabled on the float columns, producing a file that looks correct and is
    # byte-identical to one written without the split. Assert the encoding
    # landed rather than trusting the flag.
    out_path = tmp_path / "split.parquet"
    write_parquet_with_options(
        table=_table(),
        out_path=out_path,
        options=ParquetWriteOptions(byte_stream_split_floats=True),
    )
    assert "BYTE_STREAM_SPLIT" in _float_encodings(out_path)


def test_float_columns_are_not_split_by_default(tmp_path: Path):
    out_path = tmp_path / "plain.parquet"
    write_parquet_with_options(
        table=_table(),
        out_path=out_path,
        options=ParquetWriteOptions(),
    )
    assert "BYTE_STREAM_SPLIT" not in _float_encodings(out_path)


def test_values_round_trip_exactly_under_split_encoding(tmp_path: Path):
    out_path = tmp_path / "split.parquet"
    write_parquet_with_options(
        table=_table(),
        out_path=out_path,
        options=ParquetWriteOptions(
            compression="zstd", compression_level=22, byte_stream_split_floats=True
        ),
    )
    assert pl.read_parquet(out_path).equals(_frame())
