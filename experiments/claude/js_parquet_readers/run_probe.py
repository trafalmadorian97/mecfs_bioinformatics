"""
Can a browser-side JavaScript parquet reader decode BYTE_STREAM_SPLIT + zstd?

This is the load-bearing question for shipping large results tables to the docs
site as parquet instead of CSV (see docs_large_table_encoding_probe.py, which
established that BYTE_STREAM_SPLIT + zstd is what makes parquet competitive on
size). If no browser-side reader supports that encoding, the size advantage is
unreachable and CSV wins by default.

Generates parquet fixtures spanning {snappy, zstd} x {plain, BYTE_STREAM_SPLIT}
x {float64, float32} from the real PPP heritability table, then drives probe.mjs
under node to decode them with hyparquet (pure JS) and parquet-wasm
(Rust/Arrow via wasm). Decoded values are checked against a checksum of the h2
column, not merely "did it throw" --- a broken BYTE_STREAM_SPLIT implementation
would happily return plausible-looking garbage floats.

node is obtained via `pixi exec -s nodejs`, which builds a throwaway environment
and so does not touch the project's pixi.toml. npm packages are installed into a
temp directory that is removed on exit.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq

H2_PARQUET = Path(
    "assets/base_asset_store/gwas/ukbb_ppp/ppp_heritability/analysis"
    "/ppp_heritability_hapmap_3.parquet"
)

PROBE_SCRIPT = Path(__file__).parent / "probe.mjs"

# Pinned so the recorded result is attributable to specific library versions.
NPM_PACKAGES = (
    "hyparquet@1.26.2",
    "hyparquet-compressors@1.1.1",
    "parquet-wasm@0.7.2",
    "apache-arrow",
)


def write_fixture(
    frame: pl.DataFrame,
    out_path: Path,
    compression: str,
    compression_level: int | None,
    byte_stream_split: bool,
) -> None:
    """Write one parquet encoding variant.

    Dictionary encoding is disabled on exactly the float columns when
    BYTE_STREAM_SPLIT is requested, because dictionary encoding otherwise takes
    precedence and the split is silently dropped.
    """
    table = frame.to_arrow()
    float_cols = [f.name for f in table.schema if pa.types.is_floating(f.type)]
    other_cols = [f.name for f in table.schema if not pa.types.is_floating(f.type)]
    pq.write_table(
        table,
        out_path,
        compression=compression,
        compression_level=compression_level,
        use_byte_stream_split=float_cols if byte_stream_split else False,
        use_dictionary=other_cols if byte_stream_split else True,
    )


def build_fixtures(fixture_dir: Path) -> None:
    """Write the parquet fixtures and the expected-values file probe.mjs checks against."""
    assert H2_PARQUET.is_file(), f"missing input asset: {H2_PARQUET}"
    df = pl.read_parquet(H2_PARQUET)
    df32 = df.with_columns(
        pl.col(name).cast(pl.Float32)
        for name, dtype in df.schema.items()
        if dtype == pl.Float64
    )

    write_fixture(df, fixture_dir / "f64_nobss_snappy.parquet", "snappy", None, False)
    write_fixture(df, fixture_dir / "f64_nobss_zstd22.parquet", "zstd", 22, False)
    write_fixture(df, fixture_dir / "f64_bss_snappy.parquet", "snappy", None, True)
    write_fixture(df, fixture_dir / "f64_bss_zstd22.parquet", "zstd", 22, True)
    write_fixture(df32, fixture_dir / "f32_bss_zstd22.parquet", "zstd", 22, True)

    expected = {
        "n_rows": df.height,
        "cols": df.columns,
        "first": df.row(0, named=True),
        "last": df.row(-1, named=True),
        "h2_sum": float(df["h2"].sum()),
    }
    (fixture_dir / "expected.json").write_text(json.dumps(expected, indent=1))

    for path in sorted(fixture_dir.glob("*.parquet")):
        print(f"  {path.name:<28} {path.stat().st_size:>9,} bytes", flush=True)


def main() -> None:
    work_dir = Path(tempfile.mkdtemp(prefix="js_parquet_probe_"))
    fixture_dir = work_dir / "fixtures"
    fixture_dir.mkdir()

    try:
        print("fixtures:", flush=True)
        build_fixtures(fixture_dir)
        shutil.copy(PROBE_SCRIPT, work_dir / "probe.mjs")

        print("\ninstalling npm packages (throwaway env via pixi exec)...", flush=True)
        subprocess.run(
            ["pixi", "exec", "-s", "nodejs", "--", "npm", "init", "-y"],
            cwd=work_dir,
            check=True,
            capture_output=True,
        )
        subprocess.run(
            ["pixi", "exec", "-s", "nodejs", "--", "npm", "install", *NPM_PACKAGES],
            cwd=work_dir,
            check=True,
            capture_output=True,
        )

        print(flush=True)
        subprocess.run(
            ["pixi", "exec", "-s", "nodejs", "--", "node", "probe.mjs"],
            cwd=work_dir,
            check=True,
        )
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
