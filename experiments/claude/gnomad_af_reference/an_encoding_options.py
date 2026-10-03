"""
Can a different parquet encoding absorb the site-to-site jitter in gnomAD AN columns?

Throwaway follow-up to parquet_format_options.py, where the AN_g columns (Int32, near
constant per group with small jitter) cost about half the panel under the default
dictionary encoding. Only the AN columns' encoding varies here; every other column keeps
the default (dictionary where pyarrow chooses it), which was best for AF.

AN encodings compared:
- int32 dictionary (pyarrow default; the baseline);
- int32 PLAIN (no dictionary), as a control for what the codec alone does;
- float32 BYTE_STREAM_SPLIT (cast to float32 is exact: AN < 2^24);
- int32 BYTE_STREAM_SPLIT (supported for integers by recent parquet/pyarrow);
- int32 DELTA_BINARY_PACKED, parquet's encoding for slowly varying integers.

Each is written at zstd default and zstd 19. The round trip of every encoded AN column is
checked against the input.

Usage:
    pixi r python experiments/claude/gnomad_af_reference/an_encoding_options.py PANEL.parquet \
        2>&1 | tee experiments/claude/gnomad_af_reference/an_encoding_options_RELEASE.log
"""

import sys
import tempfile
import time
from pathlib import Path

import polars as pl
import pyarrow
import pyarrow.parquet
from attrs import frozen
from parquet_format_options import compressed_bytes_by_family

from mecfs_bio.build_system.task.dataframe_output import ParquetCompression


@frozen
class AnEncoding:
    name: str
    as_float32: bool
    dictionary: bool
    encoding: str | None  # explicit parquet encoding for AN columns, or None


AN_ENCODINGS = [
    AnEncoding("int32 dictionary", as_float32=False, dictionary=True, encoding=None),
    AnEncoding("int32 plain", as_float32=False, dictionary=False, encoding="PLAIN"),
    AnEncoding(
        "float32 bss", as_float32=True, dictionary=False, encoding="BYTE_STREAM_SPLIT"
    ),
    AnEncoding(
        "int32 bss", as_float32=False, dictionary=False, encoding="BYTE_STREAM_SPLIT"
    ),
    AnEncoding(
        "int32 delta",
        as_float32=False,
        dictionary=False,
        encoding="DELTA_BINARY_PACKED",
    ),
]
LEVELS: list[tuple[ParquetCompression, int | None]] = [("zstd", None), ("zstd", 19)]


def an_columns(frame: pl.DataFrame) -> list[str]:
    return [c for c in frame.columns if c.startswith("AN_")]


def write(
    frame: pl.DataFrame,
    path: Path,
    an_encoding: AnEncoding,
    compression: ParquetCompression,
    level: int | None,
) -> None:
    ans = an_columns(frame)
    if an_encoding.as_float32:
        frame = frame.with_columns(pl.col(ans).cast(pl.Float32))
    others = [c for c in frame.columns if c not in ans]
    pyarrow.parquet.write_table(
        frame.to_arrow(),
        path,
        compression=compression,
        compression_level=level,
        use_dictionary=True if an_encoding.dictionary else others,
        column_encoding=(
            {c: an_encoding.encoding for c in ans} if an_encoding.encoding else None
        ),
    )


def assert_an_round_trip(original: pl.DataFrame, path: Path) -> None:
    ans = an_columns(original)
    restored = pl.read_parquet(path, columns=ans).with_columns(
        pl.col(ans).cast(pl.Int32)
    )
    assert restored.equals(original.select(ans)), f"AN round trip failed for {path}"


def an_encodings_used(path: Path) -> str:
    metadata = pyarrow.parquet.ParquetFile(path).metadata
    row_group = metadata.row_group(0)
    for c in range(row_group.num_columns):
        chunk = row_group.column(c)
        if chunk.path_in_schema.startswith("AN_"):
            return ",".join(chunk.encodings)
    raise AssertionError("no AN column")


def main(panel_path: Path) -> None:
    frame = pl.read_parquet(panel_path)
    print(
        f"input {panel_path}: {frame.height} rows, {len(an_columns(frame))} AN columns"
    )
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for compression, level in LEVELS:
            for an_encoding in AN_ENCODINGS:
                path = Path(tmp) / f"{an_encoding.name}_{level}.parquet".replace(
                    " ", "_"
                )
                start = time.monotonic()
                write(frame, path, an_encoding, compression, level)
                write_s = time.monotonic() - start
                assert_an_round_trip(frame, path)
                families = compressed_bytes_by_family(path)
                rows.append(
                    {
                        "zstd_level": "default" if level is None else level,
                        "AN encoding": an_encoding.name,
                        "parquet encodings": an_encodings_used(path),
                        "file_mib": round(path.stat().st_size / 2**20, 2),
                        "AN_mib": round(families["AN_"] / 2**20, 2),
                        "AF_mib": round(families["AF_"] / 2**20, 2),
                        "write_s": round(write_s, 1),
                    }
                )
    with pl.Config(tbl_cols=-1, tbl_rows=-1, tbl_width_chars=250, fmt_str_lengths=60):
        print(pl.DataFrame(rows))


if __name__ == "__main__":
    assert len(sys.argv) == 2, "usage: an_encoding_options.py PANEL.parquet"
    main(Path(sys.argv[1]))
