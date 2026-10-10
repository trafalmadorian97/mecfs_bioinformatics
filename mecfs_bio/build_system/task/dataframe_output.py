"""
How a task writes a dataframe asset: the output format types, and the single
writer that interprets them.

"""

from pathlib import Path
from typing import Literal, Sequence

import narwhals
import pyarrow
import pyarrow.parquet
from attrs import frozen

from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
    DataFrameTextFormat,
)

ParquetCompression = Literal["snappy", "zstd", "gzip", "brotli", "lz4", "none"]


@frozen(slots=True)
class ParquetWriteOptions:
    """How to encode a parquet output, when the defaults are not good enough."""

    compression: ParquetCompression = "zstd"
    compression_level: int | None = None
    byte_stream_split_floats: bool = False


@frozen(slots=True)
class ParquetOutFormat:
    write_options: ParquetWriteOptions | None = None


@frozen(slots=True)
class CSVOutFormat:
    sep: str


OutFormat = ParquetOutFormat | CSVOutFormat


@frozen(slots=True)
class ParquetEncoding:
    """pyarrow parquet writer encoding arguments, as derived by parquet_encoding."""

    compression: ParquetCompression
    compression_level: int | None
    use_byte_stream_split: list[str] | bool
    use_dictionary: list[str] | bool


def parquet_encoding(
    column_names: Sequence[str],
    compression: ParquetCompression,
    compression_level: int | None,
    byte_stream_split_columns: Sequence[str],
) -> ParquetEncoding:
    """Writer encoding with BYTE_STREAM_SPLIT on the named columns and dictionary elsewhere.

    Dictionary encoding takes precedence over BYTE_STREAM_SPLIT in the parquet writer: a
    column left dictionary-enabled is written as RLE_DICTIONARY and the requested split is
    silently dropped. Dictionary encoding is therefore disabled on exactly the split
    columns and left on for the rest, where it is what makes low-cardinality columns small.
    The split applies to any fixed-width column (pyarrow 25 supports integer as well as
    floating-point types).

    Pass an empty byte_stream_split_columns to disable the split entirely.
    """
    split_columns = list(byte_stream_split_columns)
    missing = set(split_columns) - set(column_names)
    assert not missing, f"byte_stream_split_columns not in frame: {missing}"
    other_columns = [name for name in column_names if name not in set(split_columns)]
    return ParquetEncoding(
        compression=compression,
        compression_level=compression_level,
        use_byte_stream_split=split_columns if split_columns else False,
        use_dictionary=other_columns if split_columns else True,
    )


def open_parquet_writer(
    out_path: Path, schema: pyarrow.Schema, encoding: ParquetEncoding
) -> pyarrow.parquet.ParquetWriter:
    """A streaming parquet writer with the given encoding; close it (or use with)."""
    return pyarrow.parquet.ParquetWriter(
        str(out_path),
        schema,
        compression=encoding.compression,
        compression_level=encoding.compression_level,
        use_byte_stream_split=encoding.use_byte_stream_split,
        use_dictionary=encoding.use_dictionary,
    )


def write_parquet_table(
    table: pyarrow.Table,
    out_path: Path,
    compression: ParquetCompression,
    compression_level: int | None,
    byte_stream_split_columns: Sequence[str],
) -> None:
    """Write an arrow table to parquet with explicit encoding control (see parquet_encoding)."""
    encoding = parquet_encoding(
        table.schema.names,
        compression=compression,
        compression_level=compression_level,
        byte_stream_split_columns=byte_stream_split_columns,
    )
    pyarrow.parquet.write_table(
        table,
        out_path,
        compression=encoding.compression,
        compression_level=encoding.compression_level,
        use_byte_stream_split=encoding.use_byte_stream_split,
        use_dictionary=encoding.use_dictionary,
    )


def float_column_names(table: pyarrow.Table) -> list[str]:
    """Names of the floating-point columns, the ones worth byte-stream-splitting."""
    return [
        field.name for field in table.schema if pyarrow.types.is_floating(field.type)
    ]


def write_df_according_to_format(
    df: narwhals.LazyFrame, out_path: Path, out_format: OutFormat
) -> None:
    """Write a dataframe to out_path in the requested format.
    The frame's backend is preserved: narwhals dispatches to the underlying
    library's own writer
    """
    if isinstance(out_format, CSVOutFormat):
        df.collect().to_pandas().to_csv(out_path, index=False, sep=out_format.sep)
    elif isinstance(out_format, ParquetOutFormat):
        if out_format.write_options is None:
            df.sink_parquet(out_path)
        else:
            options = out_format.write_options
            table = df.collect().to_arrow()
            write_parquet_table(
                table=table,
                out_path=out_path,
                compression=options.compression,
                compression_level=options.compression_level,
                byte_stream_split_columns=(
                    float_column_names(table)
                    if options.byte_stream_split_floats
                    else []
                ),
            )
    else:
        raise ValueError(f"Unknown format {out_format}")


def get_extension_and_read_spec_from_format(
    out_format: OutFormat,
) -> tuple[str, DataFrameReadSpec]:
    if isinstance(out_format, CSVOutFormat):
        read_spec = DataFrameReadSpec(DataFrameTextFormat(separator=out_format.sep))
        if out_format.sep == "\t":
            extension = ".tsv"
        elif out_format.sep == ",":
            extension = ".csv"
        else:
            raise ValueError("Unknown sep")
    elif isinstance(out_format, ParquetOutFormat):
        read_spec = DataFrameReadSpec(DataFrameParquetFormat())
        extension = ".parquet"
    else:
        raise ValueError(f"Unknown format {out_format}")
    return extension, read_spec
