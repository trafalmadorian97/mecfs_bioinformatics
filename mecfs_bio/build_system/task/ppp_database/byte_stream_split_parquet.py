"""
Write a polars DataFrame to parquet with Zstd compression and byte-stream-split
encoding on float columns.

Byte-stream-split greatly improves compression of float mantissas but is not
exposed by polars sink_parquet/write_parquet (as of polars 1.41), so we write via
pyarrow. This is the storage format for the PPP variant index and the per-protein
beta/se files.
"""

from pathlib import Path

import polars as pl
import pyarrow.parquet as pq


def write_byte_stream_split_parquet(
    df: pl.DataFrame,
    path: Path,
    float_columns: list[str],
    compression_level: int | None = None,
) -> None:
    """Write df to path as Zstd parquet, byte-stream-split encoding float_columns.

    compression_level selects the Zstd level (None uses the pyarrow default). Pass
    an empty float_columns list to disable byte-stream-split entirely.
    """
    missing = set(float_columns) - set(df.columns)
    assert not missing, f"float_columns not in frame: {missing}"
    # Dictionary encoding takes precedence over byte-stream-split, so the split
    # columns must have dictionary disabled for BYTE_STREAM_SPLIT to be applied;
    # dictionary-encode only the remaining columns.
    dictionary_columns = [c for c in df.columns if c not in float_columns]
    pq.write_table(
        df.to_arrow(),
        str(path),
        compression="zstd",
        compression_level=compression_level,
        use_dictionary=dictionary_columns,
        use_byte_stream_split=float_columns,
    )
