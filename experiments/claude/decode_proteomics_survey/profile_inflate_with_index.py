"""
Benchmark ways to decompress one local deCODE file sequentially WHILE building a seek index
that indexed_gzip can import (GZIDX format), and check each exported index is usable.

  indexed_gzip   default / skip_crc_check + big buffers (spacing 16 MiB uncompressed)
  rapidgzip      parallelization 1, 4, 8, 16; index exported with IndexFormat.INDEXED_GZIP

Correctness check for every exported index: import it into indexed_gzip on the local file,
seek to 20 random uncompressed offsets, read 64 KiB, and compare with the same bytes from a
plain zlib decompression held as a list of 16 MiB chunks.

Usage (from outside the repo):
  pixi exec --spec python=3.12 --spec indexed_gzip --spec rapidgzip -- \
      python -I profile_inflate_with_index.py <local.txt.gz> <out_dir>
"""

import random
import sys
import time
import zlib
from pathlib import Path

import indexed_gzip
import rapidgzip

MIB = 1024 * 1024
READ = 16 * MIB


def reference_reader(path: Path):  # noqa: ANN201
    """Random access to uncompressed bytes via one zlib pass, stored as 16 MiB pieces."""
    pieces: list[bytes] = []
    d = zlib.decompressobj(31)
    buf = b""
    with open(path, "rb") as f:
        while c := f.read(4 * MIB):
            buf += d.decompress(c)
            while len(buf) >= READ:
                pieces.append(buf[:READ])
                buf = buf[READ:]
    buf += d.flush()
    if buf:
        pieces.append(buf)
    total = sum(len(p) for p in pieces)

    def get(off: int, n: int) -> bytes:
        out = b""
        while n > 0 and off < total:
            i, r = divmod(off, READ)
            take = pieces[i][r : r + n]
            out += take
            off += len(take)
            n -= len(take)
        return out

    return get, total


def check_index(path: Path, index: Path, get, total: int) -> str:  # noqa: ANN001
    rng = random.Random(0)
    with indexed_gzip.IndexedGzipFile(str(path), index_file=str(index), auto_build=False) as f:
        npoints = f.raw.npoints
        for _ in range(20):
            off = rng.randrange(0, total - 65536)
            f.seek(off)
            if f.read(65536) != get(off, 65536):
                return f"MISMATCH at {off}"
    return f"OK ({npoints} seek points)"


def run_indexed_gzip(path: Path, out: Path, **kw) -> tuple[float, int]:  # noqa: ANN003
    t = time.perf_counter()
    f = indexed_gzip.IndexedGzipFile(str(path), spacing=16 * MIB, **kw)
    n = 0
    while b := f.read(READ):
        n += len(b)
    f.raw.export_index(filename=str(out))
    f.close()
    return time.perf_counter() - t, n


def run_rapidgzip(path: Path, out: Path, parallelization: int) -> tuple[float, int]:
    t = time.perf_counter()
    n = 0
    with rapidgzip.open(str(path), parallelization=parallelization) as f:
        while b := f.read(READ):
            n += len(b)
        with open(out, "wb") as o:
            f.export_index(o, rapidgzip.IndexFormat.INDEXED_GZIP)
    return time.perf_counter() - t, n


def main() -> None:
    path, out_dir = Path(sys.argv[1]), Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)
    t = time.perf_counter()
    get, total = reference_reader(path)
    print(f"reference zlib pass: {time.perf_counter() - t:.1f}s, {total / 1e6:.0f} MB", flush=True)

    cases = [
        ("indexed_gzip default", lambda o: run_indexed_gzip(path, o)),
        ("indexed_gzip skip_crc + 16MiB readbuf", lambda o: run_indexed_gzip(
            path, o, skip_crc_check=True, readbuf_size=16 * MIB, buffer_size=64 * MIB)),
        *[(f"rapidgzip p={p}", lambda o, p=p: run_rapidgzip(path, o, p)) for p in (1, 4, 8, 16)],
    ]
    for i, (label, fn) in enumerate(cases):
        out = out_dir / f"case{i}.gzidx"
        seconds, n = fn(out)
        assert n == total, f"{label}: read {n} bytes, expected {total}"
        verdict = check_index(path, out, get, total)
        print(f"{label:40s} {seconds:6.1f}s  {total / 1e6 / seconds:5.0f} MB/s  "
              f"index {out.stat().st_size / 1e6:6.2f} MB  {verdict}", flush=True)


if __name__ == "__main__":
    main()
