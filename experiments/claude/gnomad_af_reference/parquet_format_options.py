"""
Which parquet settings minimize the size of a gnomAD allele-frequency panel?

Throwaway probe for the gnomAD AF reference brainstorm. Input is a chr21 panel written by
measure_chr21_extract.py (keys, then AF_g Float32 and AN_g Int32 per ancestry group). Each
variant rewrites it through dataframe_output.write_parquet_table and reports the file size,
the compressed bytes per column family and a full-scan read time.

Variants:
- compression codec and zstd level;
- BYTE_STREAM_SPLIT on the float (AF) columns, on or off;
- representation: AF_g Float32, or AC_g Int32 = round(AF_g * AN_g) with AF recovered as
  AC_g / AN_g at read time. The AC round trip is verified before any size is reported:
  AC must be an integer-valued count and AC / AN must reproduce the stored AF to float32
  precision.

Usage:
    pixi r python experiments/claude/gnomad_af_reference/parquet_format_options.py PANEL.parquet \
        2>&1 | tee experiments/claude/gnomad_af_reference/parquet_format_options_RELEASE.log
"""

import sys
import tempfile
import time
from pathlib import Path

import polars as pl
import pyarrow
import pyarrow.parquet
from attrs import frozen

from mecfs_bio.build_system.task.dataframe_output import (
    ParquetCompression,
    float_column_names,
    write_parquet_table,
)

# Relative AF error allowed when recovering AF as AC / AN. gnomAD prints AF with about six
# significant digits, so the recovered value differs from the printed one by rounding only.
_AF_RELATIVE_TOLERANCE = 1e-5


@frozen(slots=True)
class Variant:
    name: str
    representation: str  # "af" or "ac"
    compression: ParquetCompression
    compression_level: int | None
    byte_stream_split: bool


VARIANTS = [
    Variant("af snappy", "af", "snappy", None, False),
    Variant("af zstd(default)", "af", "zstd", None, False),
    Variant("af zstd(default) +bss", "af", "zstd", None, True),
    Variant("af zstd9", "af", "zstd", 9, False),
    Variant("af zstd9 +bss", "af", "zstd", 9, True),
    Variant("af zstd19", "af", "zstd", 19, False),
    Variant("af zstd19 +bss", "af", "zstd", 19, True),
    Variant("af brotli(default) +bss", "af", "brotli", None, True),
    Variant("ac zstd(default)", "ac", "zstd", None, False),
    Variant("ac zstd9", "ac", "zstd", 9, False),
    Variant("ac zstd19", "ac", "zstd", 19, False),
    Variant("af_only zstd(default)", "af_only", "zstd", None, False),
    Variant("af_only zstd19", "af_only", "zstd", 19, False),
    Variant("ac_deficit zstd(default)", "ac_deficit", "zstd", None, False),
    Variant("ac_deficit zstd19", "ac_deficit", "zstd", 19, False),
]


def an_deficits(frame: pl.DataFrame) -> pl.DataFrame:
    """Replace each AN_g with AN_DEF_g = max(AN_g) - AN_g, a small non-negative integer."""
    an_columns = [c for c in frame.columns if c.startswith("AN_")]
    return frame.with_columns(
        (pl.col(c).max() - pl.col(c)).alias(c.replace("AN_", "AN_DEF_"))
        for c in an_columns
    ).drop(an_columns)


def af_columns(frame: pl.DataFrame) -> list[str]:
    return [c for c in frame.columns if c.startswith("AF_")]


def to_allele_counts(frame: pl.DataFrame) -> pl.DataFrame:
    """Replace each AF_g with AC_g = round(AF_g * AN_g), after checking the round trip."""
    replacements = []
    for af in af_columns(frame):
        group = af.removeprefix("AF_")
        an = f"AN_{group}"
        # gnomAD leaves AF missing exactly when AN == 0 (no called genotypes in the group);
        # AC is then 0, and AC / AN recovers the missing AF.
        product = pl.col(af).cast(pl.Float64) * pl.col(an).cast(pl.Float64)
        replacements.append(
            product.round().fill_null(0).cast(pl.Int32).alias(f"AC_{group}")
        )
        recovered = pl.when(pl.col(an) > 0).then(pl.col(f"AC_{group}") / pl.col(an))
        check = (
            frame.with_columns(replacements[-1])
            .select(
                # AF is printed to ~6 significant digits, so AF * AN misses an integer by
                # up to ~5e-6 * AC; check it relative to the count, not absolutely.
                (
                    (product - product.round()).abs()
                    / product.round().clip(lower_bound=1)
                )
                .max()
                .alias("max_ac_rel_fraction"),
                ((recovered - pl.col(af)).abs() / pl.col(af))
                .filter(pl.col(af) > 0)
                .max()
                .alias("max_rel_err"),
                (pl.col(af).is_null() != recovered.is_null())
                .sum()
                .alias("null_mismatch"),
                pl.col(an).is_null().sum().alias("an_null"),
            )
            .row(0, named=True)
        )
        assert check["an_null"] == 0, f"{an} has nulls: {check}"
        assert check["null_mismatch"] == 0, f"{af} null except where {an} == 0: {check}"
        assert (check["max_ac_rel_fraction"] or 0.0) < _AF_RELATIVE_TOLERANCE, (
            f"{af} * {an} not integral: {check}"
        )
        assert (check["max_rel_err"] or 0.0) < _AF_RELATIVE_TOLERANCE, (
            f"AC/AN does not reproduce {af}: {check}"
        )
    return frame.with_columns(replacements).drop(af_columns(frame))


def column_family(name: str) -> str:
    if name.startswith("AN_DEF_"):
        return "AN_DEF_"
    if name.startswith(("AF_", "AC_", "AN_")):
        return name[:3]
    return name


def compressed_bytes_by_family(path: Path) -> dict[str, int]:
    metadata = pyarrow.parquet.ParquetFile(path).metadata
    totals: dict[str, int] = {}
    for rg in range(metadata.num_row_groups):
        row_group = metadata.row_group(rg)
        for c in range(row_group.num_columns):
            chunk = row_group.column(c)
            family = column_family(chunk.path_in_schema)
            totals[family] = totals.get(family, 0) + chunk.total_compressed_size
    return totals


def full_scan_seconds(path: Path) -> float:
    start = time.monotonic()
    pl.read_parquet(path)
    return time.monotonic() - start


def measure(variant: Variant, table: pyarrow.Table, out_dir: Path) -> dict[str, object]:
    path = out_dir / f"{variant.name.replace(' ', '_')}.parquet"
    write_parquet_table(
        table=table,
        out_path=path,
        compression=variant.compression,
        compression_level=variant.compression_level,
        byte_stream_split_columns=(
            float_column_names(table) if variant.byte_stream_split else []
        ),
    )
    families = compressed_bytes_by_family(path)
    return {
        "variant": variant.name,
        "file_mib": round(path.stat().st_size / 2**20, 2),
        "bytes_per_row": round(path.stat().st_size / table.num_rows, 2),
        **{f"{k}_mib": round(v / 2**20, 2) for k, v in sorted(families.items())},
        "read_s": round(min(full_scan_seconds(path) for _ in range(3)), 2),
    }


def main(panel_path: Path) -> None:
    frame = pl.read_parquet(panel_path)
    print(f"input {panel_path}: {frame.height} rows, {frame.width} columns")
    counts = to_allele_counts(frame)
    print("AC round trip verified for", len(af_columns(frame)), "groups")
    tables = {
        "af": frame.to_arrow(),
        "ac": counts.to_arrow(),
        "af_only": frame.drop(
            c for c in frame.columns if c.startswith("AN_")
        ).to_arrow(),
        "ac_deficit": an_deficits(counts).to_arrow(),
    }
    with tempfile.TemporaryDirectory() as tmp:
        rows = [measure(v, tables[v.representation], Path(tmp)) for v in VARIANTS]
    with pl.Config(tbl_cols=-1, tbl_rows=-1, tbl_width_chars=250):
        print(pl.DataFrame(rows, infer_schema_length=None))
    print(
        "NaN AF values (not nulls):",
        frame.select(
            pl.sum_horizontal(pl.col(af_columns(frame)).is_nan().sum())
        ).item(),
    )


if __name__ == "__main__":
    assert len(sys.argv) == 2, "usage: parquet_format_options.py PANEL.parquet"
    main(Path(sys.argv[1]))
