from typing import Literal, Mapping, Sequence

import polars as pl
from attrs import field, frozen

Compression = Literal["gzip"]


@frozen
class DataFrameTextFormat:
    separator: str
    null_values: Sequence[str] | None = None
    schema_overrides: Mapping[str, pl.DataType] = field(factory=dict)
    column_names: list[str] | None = None
    has_header: bool = True
    skip_rows: int = 0


@frozen
class DataFrameParquetFormat:
    pass


@frozen
class DataFrameWhiteSpaceSepTextFormat:
    comment_code: str


DataFrameFormat = (
    DataFrameTextFormat | DataFrameParquetFormat | DataFrameWhiteSpaceSepTextFormat
)


@frozen
class DataFrameReadSpec:
    """
    Specifies how a file containing a dataframe should be read.
    Allows client code to operate without concern for the specifics of the dataframe format.
    """

    format: DataFrameFormat
