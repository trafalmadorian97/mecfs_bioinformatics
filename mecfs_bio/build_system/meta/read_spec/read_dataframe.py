from pathlib import Path
from typing import Literal

import narwhals as nw
import pandas as pd
import polars as pl

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.base_meta import FileMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
    DataFrameTextFormat,
    DataFrameWhiteSpaceSepTextFormat,
)

ValidBackend = Literal["ibis", "polars", "duckdb"]
import structlog

logger = structlog.get_logger()


def scan_dataframe(
    path: Path, spec: DataFrameReadSpec, parquet_backend: ValidBackend = "polars"
) -> nw.LazyFrame:
    if isinstance(spec.format, DataFrameParquetFormat):
        logger.debug(f"Scanning with backend {parquet_backend}")
        return nw.scan_parquet(path, backend=parquet_backend)
    if isinstance(spec.format, DataFrameTextFormat):
        if spec.format.column_names is not None:
            col_list: list[str] = spec.format.column_names
            col_func = lambda x: col_list
        else:
            col_func = None
        if parquet_backend == "polars":
            polars_scan = pl.scan_csv(
                path,
                separator=spec.format.separator,
                null_values=spec.format.null_values,
                schema_overrides=spec.format.schema_overrides,
                with_column_names=col_func,
                has_header=spec.format.has_header,
                skip_rows=spec.format.skip_rows,
                comment_prefix=spec.format.comment_char,
            )
            return nw.from_native(polars_scan)
        raise ValueError("Only polars backend can be used to read text files")
    if isinstance(spec.format, DataFrameWhiteSpaceSepTextFormat):
        return nw.from_native(
            pl.from_pandas(
                pd.read_csv(path, sep=r"\s+", comment=spec.format.comment_code)
            )
        ).lazy()

    raise ValueError("Unknown format")


def _scan_dataframe_asset(
    asset: FileAsset, meta: FileMeta, parquet_backend: ValidBackend
) -> nw.LazyFrame:
    read_spec = meta.read_spec()
    assert read_spec is not None
    return scan_dataframe(asset.path, read_spec, parquet_backend)


def scan_dataframe_asset(
    asset: Asset, meta: Meta, parquet_backend: ValidBackend = "polars"
) -> nw.LazyFrame:
    """
    Use the information in an Asset's metadata to read it as a DataFrame
    """
    assert isinstance(asset, FileAsset)
    assert isinstance(meta, FileMeta)
    return _scan_dataframe_asset(asset, meta, parquet_backend)


