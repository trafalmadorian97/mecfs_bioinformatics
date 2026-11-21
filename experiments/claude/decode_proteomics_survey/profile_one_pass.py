"""
Profile the pieces of the mockup's 51 s single decompression pass over one local deCODE
SomaScan file (see mockup_single_protein.py), each timed in isolation:

  inflate_indexed_gzip  indexed_gzip sequential read, building a 16 MiB-spaced seek index
  inflate_zlib          stdlib zlib.decompressobj
  inflate_isal          python-isal (ISA-L) igzip decompressobj
  newline_scan          numpy: positions of every newline (row boundaries)
  polars_full_parse     polars read_csv of every line (what the mockup does)
  hapmap3_join          the semi-join to HapMap3 on already-parsed blocks
  gather_parse          parse ONLY the lines whose row numbers the full parse found in HapMap3:
                        numpy newline scan -> slice those lines -> polars parse; checked equal
                        to the full-parse result

Decompression for the parse stages comes from zlib and is excluded from their timings.

Usage (from outside the repo):
  pixi exec --spec python=3.12 --spec indexed_gzip --spec polars --spec python-isal \
      --spec numpy -- python -I profile_one_pass.py <local.txt.gz> <hapmap3.parquet>
"""

import io
import sys
import time
import zlib
from collections.abc import Iterator
from pathlib import Path

import indexed_gzip
import numpy as np
import polars as pl
from isal import isal_zlib

MIB = 1024 * 1024
BLOCK = 16 * MIB
HEADER = [
    "Chrom", "Pos", "Name", "rsids", "effectAllele", "otherAllele", "Beta", "Pval",
    "minus_log10_pval", "SE", "N", "ImpMAF",
]
TYPES = [
    pl.String, pl.Int64, pl.String, pl.String, pl.String, pl.String, pl.Float64,
    pl.Float64, pl.Float64, pl.Float64, pl.Float64, pl.Float64,
]
SCHEMA = pl.Schema(dict(zip(HEADER, TYPES, strict=True)))
KEEP = ["Chrom", "Pos", "Name", "effectAllele", "otherAllele", "Beta", "SE", "N", "ImpMAF"]


class Timer:
    def __init__(self) -> None:
        self.totals: dict[str, float] = {}

    def add(self, name: str, seconds: float) -> None:
        self.totals[name] = self.totals.get(name, 0.0) + seconds


def inflate_indexed_gzip(path: Path) -> Iterator[bytes]:
    f = indexed_gzip.IndexedGzipFile(str(path), spacing=16 * MIB, readbuf_size=4 * MIB)
    while b := f.read(BLOCK):
        yield b
    f.close()


def inflate_with(path: Path, decompressobj) -> Iterator[bytes]:  # noqa: ANN001
    d = decompressobj
    with open(path, "rb") as f:
        while c := f.read(4 * MIB):
            yield d.decompress(c)
        yield d.flush()


def line_blocks(chunks: Iterator[bytes]) -> Iterator[bytes]:
    """Re-cut decompressed chunks into blocks of whole lines; drops the header line."""
    pending = b""
    header_done = False
    for c in chunks:
        pending += c
        if len(pending) < BLOCK:
            continue
        cut = pending.rindex(b"\n") + 1
        block, pending = pending[:cut], pending[cut:]
        if not header_done:
            block = block[block.index(b"\n") + 1 :]
            header_done = True
        yield block
    if pending:
        yield pending if pending.endswith(b"\n") else pending + b"\n"


def parse(text: bytes) -> pl.DataFrame:
    return pl.read_csv(
        io.BytesIO(text), separator="\t", has_header=False, schema=SCHEMA, null_values=["."]
    ).select(KEEP)


def pair(a: str, b: str) -> list[pl.Expr]:
    lo = pl.when(pl.col(a) <= pl.col(b)).then(pl.col(a)).otherwise(pl.col(b))
    hi = pl.when(pl.col(a) <= pl.col(b)).then(pl.col(b)).otherwise(pl.col(a))
    return [lo.alias("_a"), hi.alias("_b")]


def hm3_join(rows: pl.DataFrame, hm3: pl.DataFrame) -> pl.DataFrame:
    keyed = rows.with_columns(
        pl.col("Chrom").str.strip_prefix("chr").replace({"X": "23"}).cast(pl.Int32).alias("CHR"),
        pl.col("Pos").cast(pl.Int32).alias("POS"),
        *pair("effectAllele", "otherAllele"),
    )
    return keyed.join(hm3, on=["CHR", "POS", "_a", "_b"], how="semi").drop("CHR", "POS", "_a", "_b")


def time_inflate(name: str, chunks: Iterator[bytes], timer: Timer) -> int:
    t = time.perf_counter()
    n = sum(len(c) for c in chunks)
    timer.add(name, time.perf_counter() - t)
    return n


def main() -> None:
    path, hm3_path = Path(sys.argv[1]), Path(sys.argv[2])
    hm3 = pl.read_parquet(hm3_path).select("CHR", "POS", *pair("EA", "NEA"))
    timer = Timer()

    for name, chunks in [
        ("inflate_indexed_gzip", inflate_indexed_gzip(path)),
        ("inflate_zlib", inflate_with(path, zlib.decompressobj(31))),
        ("inflate_isal", inflate_with(path, isal_zlib.decompressobj(31))),
    ]:
        n = time_inflate(name, chunks, timer)
        print(f"{name}: {timer.totals[name]:.1f}s for {n / 1e6:.0f} MB "
              f"({n / 1e6 / timer.totals[name]:.0f} MB/s uncompressed)", flush=True)

    # Full parse + join, recording the global row numbers of the HapMap3 rows.
    row_base = 0
    kept_rows: list[np.ndarray] = []
    full_parts: list[pl.DataFrame] = []
    blocks = 0
    for block in line_blocks(inflate_with(path, isal_zlib.decompressobj(31))):
        blocks += 1
        t = time.perf_counter()
        arr = np.frombuffer(block, dtype=np.uint8)
        newlines = np.flatnonzero(arr == 10)
        timer.add("newline_scan", time.perf_counter() - t)

        t = time.perf_counter()
        rows = parse(block)
        timer.add("polars_full_parse", time.perf_counter() - t)
        assert rows.height == len(newlines), (rows.height, len(newlines))

        t = time.perf_counter()
        kept = hm3_join(rows.with_row_index("_r"), hm3)
        timer.add("hapmap3_join", time.perf_counter() - t)
        kept_rows.append(kept["_r"].to_numpy().astype(np.int64) + row_base)
        full_parts.append(kept.drop("_r"))
        row_base += rows.height
    full = pl.concat(full_parts)
    targets = np.concatenate(kept_rows)
    print(f"full parse: {row_base:,} rows in {blocks} blocks; {full.height:,} HapMap3 rows", flush=True)

    # Gather: parse only the target lines, located by row number via the newline scan.
    row_base = 0
    gather_parts: list[pl.DataFrame] = []
    for block in line_blocks(inflate_with(path, isal_zlib.decompressobj(31))):
        t = time.perf_counter()
        arr = np.frombuffer(block, dtype=np.uint8)
        ends = np.flatnonzero(arr == 10)
        starts = np.concatenate(([0], ends[:-1] + 1))
        lo = np.searchsorted(targets, row_base)
        hi = np.searchsorted(targets, row_base + len(ends))
        local = targets[lo:hi] - row_base
        selected = b"".join(block[s : e + 1] for s, e in zip(starts[local], ends[local], strict=True))
        part = parse(selected)
        timer.add("gather_parse", time.perf_counter() - t)
        gather_parts.append(part)
        row_base += len(ends)
    gathered = pl.concat(gather_parts)
    print(f"gather parse equals full parse+join: {gathered.equals(full)}", flush=True)

    print("\nTIMINGS (s)")
    for k, v in timer.totals.items():
        print(f"  {k:22s} {v:6.1f}")


if __name__ == "__main__":
    main()
