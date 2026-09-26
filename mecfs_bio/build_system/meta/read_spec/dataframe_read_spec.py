from typing import Literal, Mapping, Sequence

import polars as pl
from attrs import field, frozen

Compression = Literal["gzip"]


@frozen(slots=True)
class DataFrameTextFormat:
    separator: str
    null_values: Sequence[str] | None = None
    schema_overrides: Mapping[str, pl.DataType] = field(factory=dict)
    column_names: list[str] | None = None
    has_header: bool = True
    skip_rows: int = 0
    comment_char: str | None = None


@frozen(slots=True)
class DataFrameParquetFormat:
    pass


@frozen(slots=True)
class DataFrameWhiteSpaceSepTextFormat:
    comment_code: str
    col_names: list[str] | None = None


DataFrameFormat = (
    DataFrameTextFormat | DataFrameParquetFormat | DataFrameWhiteSpaceSepTextFormat
)


@frozen(slots=True)
class DataFrameReadSpec:
    """
    Specifies how a file containing a dataframe should be read.
    Allows client code to operate without concern for the specifics of the dataframe format.
    """

    format: DataFrameFormat
