from pathlib import Path

import narwhals
import pandas as pd
import polars as pl
import pyarrow.parquet as pq

from mecfs_bio.build_system.task.dataframe_output import (
    CSVOutFormat,
    ParquetOutFormat,
    ParquetWriteOptions,
    write_df_according_to_format,
    write_parquet_table,
)


def _frame() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "label": ["a", "b", "a", "b"],
            "count": [1, 2, 3, 4],
            "value": [0.125, 0.25, 0.5, 1.0],
        }
    )


def _lazy(backend: str) -> narwhals.LazyFrame:
    frame = _frame()
    native = frame.to_pandas() if backend == "pandas" else frame
    return narwhals.from_native(native).lazy()


def _encodings(path: Path, column: str) -> set[str]:
    row_group = pq.ParquetFile(path).metadata.row_group(0)
    columns = [row_group.column(i) for i in range(row_group.num_columns)]
    return set(next(c for c in columns if c.path_in_schema == column).encodings)


def test_byte_stream_split_is_applied_to_named_columns(tmp_path: Path):
    # Dictionary encoding silently wins over BYTE_STREAM_SPLIT if it is left
    # enabled on the split columns, producing a file that looks correct and is
    # byte-identical to one written without the split. Assert the encoding
    # landed rather than trusting the flag.
    out_path = tmp_path / "split.parquet"
    write_parquet_table(
        table=_frame().to_arrow(),
        out_path=out_path,
        compression="zstd",
        compression_level=None,
        byte_stream_split_columns=["value"],
    )
    assert "BYTE_STREAM_SPLIT" in _encodings(out_path, "value")


def test_empty_split_column_list_disables_the_split(tmp_path: Path):
    out_path = tmp_path / "plain.parquet"
    write_parquet_table(
        table=_frame().to_arrow(),
        out_path=out_path,
        compression="zstd",
        compression_level=None,
        byte_stream_split_columns=[],
    )
    assert "BYTE_STREAM_SPLIT" not in _encodings(out_path, "value")


def test_write_options_split_floats_only(tmp_path: Path):
    out_path = tmp_path / "opts.parquet"
    write_df_according_to_format(
        df=_lazy("polars"),
        out_path=out_path,
        out_format=ParquetOutFormat(
            ParquetWriteOptions(
                compression="zstd", compression_level=22, byte_stream_split_floats=True
            )
        ),
    )
    assert "BYTE_STREAM_SPLIT" in _encodings(out_path, "value")
    # Non-float columns keep dictionary encoding, which is what makes
    # low-cardinality string columns small.
    assert "RLE_DICTIONARY" in _encodings(out_path, "label")
    assert pl.read_parquet(out_path).equals(_frame())


def test_parquet_without_options_round_trips(tmp_path: Path):
    out_path = tmp_path / "default.parquet"
    write_df_according_to_format(
        df=_lazy("polars"), out_path=out_path, out_format=ParquetOutFormat()
    )
    assert pl.read_parquet(out_path).equals(_frame())


def test_csv_format_uses_the_requested_separator(tmp_path: Path):
    out_path = tmp_path / "out.tsv"
    write_df_according_to_format(
        df=_lazy("polars"), out_path=out_path, out_format=CSVOutFormat(sep="\t")
    )
    assert out_path.read_text().splitlines()[0] == "label\tcount\tvalue"


def test_pandas_backed_parquet_matches_pandas_own_writer(tmp_path: Path):
    # narwhals dispatches to the backend's own writer, so routing a pandas-backed
    # task through the shared helper must not change a single byte of its output.
    via_helper = tmp_path / "helper.parquet"
    via_pandas = tmp_path / "pandas.parquet"
    write_df_according_to_format(
        df=_lazy("pandas"), out_path=via_helper, out_format=ParquetOutFormat()
    )
    _frame().to_pandas().to_parquet(via_pandas)
    assert via_helper.read_bytes() == via_pandas.read_bytes()


def test_pandas_backed_csv_matches_pandas_own_writer(tmp_path: Path):
    via_helper = tmp_path / "helper.csv"
    via_pandas = tmp_path / "pandas.csv"
    write_df_according_to_format(
        df=_lazy("pandas"), out_path=via_helper, out_format=CSVOutFormat(sep=",")
    )
    _frame().to_pandas().to_csv(via_pandas, index=False, sep=",")
    assert via_helper.read_text() == via_pandas.read_text()


def test_pandas_backed_frame_with_gap_index_preserves_pandas_behaviour(tmp_path: Path):
    # A filtered pandas frame carries a non-contiguous index, which pandas
    # persists as __index_level_0__. The helper must not silently drop it.
    filtered = _frame().to_pandas().query("count > 2")
    via_helper = tmp_path / "helper.parquet"
    via_pandas = tmp_path / "pandas.parquet"
    write_df_according_to_format(
        df=narwhals.from_native(filtered).lazy(),
        out_path=via_helper,
        out_format=ParquetOutFormat(),
    )
    filtered.to_parquet(via_pandas)
    assert via_helper.read_bytes() == via_pandas.read_bytes()


def test_pandas_frame_is_accepted_via_from_native(tmp_path: Path):
    out_path = tmp_path / "out.parquet"
    write_df_according_to_format(
        df=narwhals.from_native(pd.DataFrame({"x": [1, 2]})).lazy(),
        out_path=out_path,
        out_format=ParquetOutFormat(),
    )
    assert pl.read_parquet(out_path)["x"].to_list() == [1, 2]
