"""
How small can the PPP per-protein heritability table get on disk?

Context: we are deciding how to ship a ~5.9k-row results table to the docs site
(see the docs figure system). A CSV is trivially consumable by a JavaScript
table widget but is large; a parquet is smaller but needs a JS parquet reader.
Before paying for that dependency we want to know how much parquet actually
buys once it is encoded properly, rather than at polars' defaults.

Sweep, on the real h2 table:
  - container: csv / csv.gz / parquet
  - float width: float64 vs float32
  - rounding: full precision vs rounded to 5 decimal places
  - compression: uncompressed / snappy / zstd (levels 3 and 22)
  - encoding: default (plain/dictionary) vs BYTE_STREAM_SPLIT on float columns

BYTE_STREAM_SPLIT is the interesting one: it transposes the mantissa/exponent
bytes of a float column so that like-magnitude bytes sit together, which
usually makes a general-purpose compressor far more effective on float data.

Writes every candidate to a temp dir, measures the on-disk size, and prints a
sorted table. Also verifies that each parquet variant round-trips to the same
values (within the tolerance implied by its float width / rounding) so we are
not comparing a lossy encoding against a lossless one without noticing.
"""

from __future__ import annotations

import gzip
import shutil
import sys
import tempfile
from pathlib import Path

import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq

# The analysis asset produced by the PPP heritability pipeline.
H2_PARQUET = Path(
    "assets/base_asset_store/gwas/ukbb_ppp/ppp_heritability/analysis"
    "/ppp_heritability_hapmap_3.parquet"
)

# Decimal places used when rounding. The h2 standard errors are ~1e-2, so five
# decimals is already well past the point of scientific meaning; this matches
# the precision the committed markdown table displays.
ROUND_DP = 5

# zstd levels worth contrasting: 3 is the usual default, 22 is the maximum and
# tells us how much headroom brute-force compression has.
ZSTD_LEVELS = (3, 22)


def float_columns(df: pl.DataFrame) -> list[str]:
    """Names of the floating-point columns, which are the ones worth re-encoding."""
    return [name for name, dtype in df.schema.items() if dtype in (pl.Float64, pl.Float32)]


def round_floats(df: pl.DataFrame, decimals: int) -> pl.DataFrame:
    """Round every float column to a fixed number of decimal places."""
    return df.with_columns(pl.col(name).round(decimals) for name in float_columns(df))


def to_float32(df: pl.DataFrame) -> pl.DataFrame:
    """Narrow every float column to 32 bits."""
    return df.with_columns(pl.col(name).cast(pl.Float32) for name in float_columns(df))


def write_parquet_variant(
    table: pa.Table,
    out_path: Path,
    compression: str,
    compression_level: int | None,
    byte_stream_split: bool,
) -> None:
    """Write one parquet encoding variant.

    byte_stream_split is applied only to the float columns; applying it to the
    string or integer columns is either unsupported or counterproductive, since
    those benefit from dictionary and RLE encoding instead.

    Critically, dictionary encoding takes precedence over BYTE_STREAM_SPLIT in
    the parquet writer: a column left dictionary-enabled comes out
    RLE_DICTIONARY and the split is silently ignored. So dictionary encoding is
    disabled *per column* on exactly the float columns, and left on for the
    low-cardinality string columns (variant_set has 2 distinct values, and
    gene/oid each repeat across the two variant sets), which is where the
    string-side savings come from.
    """
    float_cols = [
        field.name
        for field in table.schema
        if pa.types.is_floating(field.type)
    ]
    non_float_cols = [
        field.name
        for field in table.schema
        if not pa.types.is_floating(field.type)
    ]
    pq.write_table(
        table,
        out_path,
        compression=compression,
        compression_level=compression_level,
        use_byte_stream_split=float_cols if byte_stream_split else False,
        use_dictionary=non_float_cols if byte_stream_split else True,
    )


def verify_roundtrip(out_path: Path, expected: pl.DataFrame) -> None:
    """Assert the written parquet reads back to exactly the frame we handed it.

    Every lossy step (rounding, float32 narrowing) is applied to the frame
    *before* it is written, so the write itself must be lossless regardless of
    which compression or encoding variant produced it.
    """
    actual = pl.read_parquet(out_path)
    assert actual.equals(expected), f"round-trip mismatch for {out_path.name}"


def verify_encoding(out_path: Path, byte_stream_split: bool) -> None:
    """Assert the float columns really got the encoding we asked for.

    Guards against the failure mode where dictionary encoding silently wins and
    BYTE_STREAM_SPLIT is dropped, which makes the split and non-split variants
    come out byte-identical and the whole comparison meaningless.
    """
    metadata = pq.ParquetFile(out_path).metadata
    schema = pq.ParquetFile(out_path).schema_arrow
    float_cols = {field.name for field in schema if pa.types.is_floating(field.type)}

    row_group = metadata.row_group(0)
    for i in range(row_group.num_columns):
        column = row_group.column(i)
        if column.path_in_schema not in float_cols:
            continue
        encodings = set(column.encodings)
        if byte_stream_split:
            assert "BYTE_STREAM_SPLIT" in encodings, (
                f"{out_path.name}: expected BYTE_STREAM_SPLIT on "
                f"'{column.path_in_schema}', got {sorted(encodings)}"
            )
        else:
            assert "BYTE_STREAM_SPLIT" not in encodings, (
                f"{out_path.name}: unexpected BYTE_STREAM_SPLIT on "
                f"'{column.path_in_schema}'"
            )


def main() -> None:
    assert H2_PARQUET.is_file(), f"missing input asset: {H2_PARQUET}"

    base = pl.read_parquet(H2_PARQUET)
    print(f"input: {H2_PARQUET}")
    print(f"shape: {base.shape}")
    print(f"float columns: {float_columns(base)}")
    print()

    # The four (precision, width) combinations we care about, cheapest last.
    frames: dict[str, pl.DataFrame] = {
        "f64/full": base,
        "f64/round5": round_floats(base, ROUND_DP),
        "f32/full": to_float32(base),
        "f32/round5": to_float32(round_floats(base, ROUND_DP)),
    }

    results: list[tuple[str, int]] = []
    work_dir = Path(tempfile.mkdtemp(prefix="h2_encoding_probe_"))

    try:
        for frame_label, frame in frames.items():
            # --- CSV baselines (only meaningful for the float64 frames, but we
            # measure all of them since CSV serialises the rounded/narrowed
            # values differently). ---
            csv_path = work_dir / f"{frame_label.replace('/', '_')}.csv"
            frame.write_csv(csv_path)
            results.append((f"csv            [{frame_label}]", csv_path.stat().st_size))

            gz_path = csv_path.with_suffix(".csv.gz")
            gz_path.write_bytes(gzip.compress(csv_path.read_bytes(), compresslevel=9))
            results.append((f"csv.gz         [{frame_label}]", gz_path.stat().st_size))

            # --- Parquet variants. ---
            arrow_table = frame.to_arrow()

            variants: list[tuple[str, str, int | None]] = [
                ("uncompressed", "none", None),
                ("snappy", "snappy", None),
                *[(f"zstd-{lvl}", "zstd", lvl) for lvl in ZSTD_LEVELS],
            ]

            for variant_label, compression, level in variants:
                for split in (False, True):
                    suffix = "+bss" if split else ""
                    name = f"{frame_label.replace('/', '_')}_{variant_label}{suffix}"
                    path = work_dir / f"{name}.parquet"
                    write_parquet_variant(
                        arrow_table,
                        path,
                        compression=compression,
                        compression_level=level,
                        byte_stream_split=split,
                    )
                    verify_roundtrip(path, frame)
                    verify_encoding(path, byte_stream_split=split)
                    label = f"parquet {variant_label}{suffix}".ljust(15)
                    results.append((f"{label}[{frame_label}]", path.stat().st_size))

        # --- Report, smallest first. ---
        results.sort(key=lambda row: row[1])
        smallest = results[0][1]
        csv_full = next(size for label, size in results if label.startswith("csv  ") and "f64/full" in label)

        print(f"{'variant':<34} {'bytes':>9} {'KB':>8} {'vs csv f64':>11} {'vs best':>8}")
        print("-" * 74)
        for label, size in results:
            print(
                f"{label:<34} {size:>9,} {size / 1024:>7.0f}K "
                f"{csv_full / size:>10.1f}x {size / smallest:>7.2f}x"
            )
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
