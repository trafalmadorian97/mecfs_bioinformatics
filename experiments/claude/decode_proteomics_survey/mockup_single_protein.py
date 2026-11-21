"""
End-to-end mockup of the planned deCODE SomaScan proteomics database build and coloc window
reader, rewired after profiling (profile_*.py): rapidgzip builds the seek index, and the
HapMap3 rows are gathered by row number from a shared variant index instead of parsing and
joining all 33M rows.

Two modes, run as separate processes so the protein mode's peak memory is its own:

  template  (once per SomaScan group) download a template protein file, fully parse it, keep
            the HapMap3 rows, and write the shared variant index: row_number, Name, Chrom,
            Pos, effectAllele, otherAllele, plus the file's total row count.

  protein   the per-protein build, then the coloc window reader:
     1. presign a URL with boto3; read deCODE's md5 from the object metadata (HEAD)
     2. download with many parallel HTTP Range requests into a preallocated file; check md5
     3. ONE rapidgzip pass (parallel decompression) that
          - builds the seek index (one seek point per --chunk-size-mib of compressed input),
          - records the position map (first line at or after each 1 MiB of uncompressed text),
          - locates row boundaries with a numpy newline scan and parses ONLY the shared
            index's rows, checking every gathered Name against the index and the file's total
            row count against the template's
     4. write the slim file (Beta, SE as float32, in index row order; zstd + byte-stream-
        split), the zstd-compressed GZIDX seek index (rapidgzip export, indexed_gzip format)
        and the position map
     5. remote window reads (+-500 kb around the strongest index variant, and one late-file
        chr9 window) through indexed_gzip over a prefetched remote file-like object
     6. each window checked against a brute-force gzip -dc | awk scan of the local file

Needs AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY for the deCODE bucket in the environment.
Run from a directory outside the repo:

  ENV="--spec python=3.12 --spec indexed_gzip --spec rapidgzip --spec polars --spec numpy
       --spec pyarrow --spec zstandard --spec boto3"
  pixi exec $ENV -- python -I mockup_single_protein.py template --workdir W --hapmap3 H
  pixi exec $ENV -- python -I mockup_single_protein.py protein --workdir W
"""

import argparse
import base64
import hashlib
import io
import os
import resource
import subprocess
import sys
import threading
import time
import urllib.request
from bisect import bisect_right
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import boto3
import indexed_gzip
import numpy as np
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq
import rapidgzip
import zstandard
from botocore.config import Config

ENDPOINT = "https://s3-ext.decode.is"
BUCKET = "largescaleplasma-2023"
TEMPLATE_KEY = "final_somascan_smp/Proteomics_SMP_PC0_10000_28_CRYBB2_CRBB2_10032022.txt.gz"
PROTEIN_KEY = "final_somascan_smp/Proteomics_SMP_PC0_10076_1_AP4M1_AP4M1_10032022.txt.gz"
MIB = 1024 * 1024
READ = 16 * MIB
POSMAP_STEP = 1 * MIB
HALF_WINDOW_BP = 500_000
HEADER = [
    "Chrom", "Pos", "Name", "rsids", "effectAllele", "otherAllele", "Beta", "Pval",
    "minus_log10_pval", "SE", "N", "ImpMAF",
]
TYPES = [
    pl.String, pl.Int64, pl.String, pl.String, pl.String, pl.String, pl.Float64,
    pl.Float64, pl.Float64, pl.Float64, pl.Float64, pl.Float64,
]
FULL_SCHEMA = pl.Schema(dict(zip(HEADER, TYPES, strict=True)))
KEEP = ["Chrom", "Pos", "Name", "effectAllele", "otherAllele", "Beta", "SE", "N", "ImpMAF"]
INDEX_COLS = ["row_number", "Name", "Chrom", "Pos", "effectAllele", "otherAllele"]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def peak_rss_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def chrom_code(chrom: str) -> int:
    c = chrom.removeprefix("chr")
    return 23 if c == "X" else int(c)


def chrom_code_expr(col: str) -> pl.Expr:
    return pl.col(col).str.strip_prefix("chr").replace({"X": "23"}).cast(pl.Int32)


def parse_lines(text: bytes) -> pl.DataFrame:
    """Parse header-less lines of a deCODE sumstats file into the KEEP columns."""
    if not text:
        return pl.DataFrame(schema=FULL_SCHEMA).select(KEEP)
    return pl.read_csv(
        io.BytesIO(text), separator="\t", has_header=False, schema=FULL_SCHEMA,
        null_values=["."],
    ).select(KEEP)


# ---------------------------------------------------------------- S3 + range requests


@dataclass(frozen=True)
class RemoteObject:
    url: str
    size: int
    md5_hex: str


def remote_object(key: str) -> RemoteObject:
    s3 = boto3.client(
        "s3", endpoint_url=ENDPOINT, region_name="us-east-1",
        config=Config(signature_version="s3v4", s3={"addressing_style": "path"}),
    )
    head = s3.head_object(Bucket=BUCKET, Key=key)
    md5_hex = base64.b64decode(head["Metadata"]["md5chksum"]).hex()
    url = s3.generate_presigned_url(
        "get_object", Params={"Bucket": BUCKET, "Key": key}, ExpiresIn=6 * 3600
    )
    return RemoteObject(url=url, size=int(head["ContentLength"]), md5_hex=md5_hex)


@dataclass
class FetchStats:
    requests: int = 0
    bytes: int = 0
    lock: threading.Lock = field(default_factory=threading.Lock)

    def add(self, n: int) -> None:
        with self.lock:
            self.requests += 1
            self.bytes += n


def fetch_range(url: str, start: int, end: int, stats: FetchStats) -> bytes:
    """Bytes [start, end] (inclusive) of url, with retries."""
    for attempt in range(8):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(req, timeout=120) as r:
                assert r.status == 206, f"expected 206, got {r.status}"
                data = r.read()
            assert len(data) == end - start + 1, f"short read at {start}"
            stats.add(len(data))
            return data
        except (OSError, AssertionError) as e:
            print(f"  retry {attempt} for range {start}: {e}", file=sys.stderr)
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"range {start}-{end} failed")


def chunk_bounds(start: int, stop: int, chunk: int) -> list[tuple[int, int]]:
    """Inclusive (start, end) pieces covering [start, stop)."""
    return [(s, min(s + chunk, stop) - 1) for s in range(start, stop, chunk)]


def parallel_download(obj: RemoteObject, out: Path, threads: int, chunk: int) -> FetchStats:
    stats = FetchStats()
    fd = os.open(out, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    try:
        os.ftruncate(fd, obj.size)

        def one(bounds: tuple[int, int]) -> None:
            os.pwrite(fd, fetch_range(obj.url, bounds[0], bounds[1], stats), bounds[0])

        with ThreadPoolExecutor(threads) as pool:
            list(pool.map(one, chunk_bounds(0, obj.size, chunk)))
    finally:
        os.close(fd)
    return stats


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        while b := f.read(16 * MIB):
            h.update(b)
    return h.hexdigest()


def download_verified(key: str, local: Path, threads: int, chunk_mib: int,
                      timings: dict[str, float]) -> RemoteObject:
    t = time.time()
    obj = remote_object(key)
    timings["presign+head"] = time.time() - t
    log(f"  {key}: {obj.size / 1e6:.0f} MB, deCODE md5 {obj.md5_hex}")
    t = time.time()
    dl = parallel_download(obj, local, threads=threads, chunk=chunk_mib * MIB)
    timings["download"] = time.time() - t
    log(f"  downloaded {dl.bytes / 1e6:.0f} MB in {dl.requests} requests, "
        f"{timings['download']:.0f}s = {dl.bytes / 1e6 / timings['download']:.1f} MB/s")
    t = time.time()
    got = md5_of(local)
    timings["md5"] = time.time() - t
    assert got == obj.md5_hex, f"md5 mismatch: {got} != {obj.md5_hex}"
    log("  md5 OK")
    return obj


# ---------------------------------------------------------------- line blocks


@dataclass(frozen=True)
class LineBlock:
    """Whole lines of the file: text, its uncompressed offset, and its first row number."""

    text: bytes
    offset: int
    first_row: int
    line_ends: np.ndarray  # positions of the newline ending each line, relative to text


def line_blocks(f) -> Iterator[LineBlock]:  # noqa: ANN001
    """Cut a sequential decompressed stream into blocks of whole data lines (header checked
    and skipped), locating every line end with a numpy newline scan."""
    first = f.read(READ)
    header_end = first.index(b"\n") + 1
    assert first[: header_end - 1].decode().split("\t") == HEADER, "unexpected header"
    pending, offset, row = first[header_end:], header_end, 0
    eof = False
    while not eof:
        more = f.read(READ)
        eof = not more
        pending += more
        if not eof and len(pending) < READ:
            continue
        if eof and pending and not pending.endswith(b"\n"):
            pending += b"\n"
        cut = len(pending) if eof else pending.rindex(b"\n") + 1
        text, pending = pending[:cut], pending[cut:]
        if text:
            ends = np.flatnonzero(np.frombuffer(text, dtype=np.uint8) == 10)
            yield LineBlock(text=text, offset=offset, first_row=row, line_ends=ends)
            offset += len(text)
            row += len(ends)


def posmap_rows(block: LineBlock, next_mark: int) -> tuple[list[tuple[int, int, int, int]], int]:
    """Position-map rows for the first line starting at or after each POSMAP_STEP mark."""
    starts = np.concatenate(([0], block.line_ends[:-1] + 1))
    rows: list[tuple[int, int, int, int]] = []
    while next_mark < block.offset + len(block.text):
        i = int(np.searchsorted(starts, next_mark - block.offset, side="left"))
        if i >= len(starts):
            break
        s = int(starts[i])
        fields = block.text[s : s + 64].split(b"\t", 2)
        rows.append((chrom_code(fields[0].decode()), int(fields[1]), block.offset + s,
                     block.first_row + i))
        next_mark = block.offset + s + POSMAP_STEP
    return rows, next_mark


def gather_lines(block: LineBlock, rows: np.ndarray) -> bytes:
    """The text of the given (global) row numbers, all inside this block."""
    local = rows - block.first_row
    ends = block.line_ends[local]
    starts = np.where(local == 0, 0, block.line_ends[np.maximum(local - 1, 0)] + 1)
    return b"".join(block.text[s : e + 1] for s, e in zip(starts.tolist(), ends.tolist(), strict=True))


# ---------------------------------------------------------------- template mode


def sorted_allele_pair(a: str, b: str) -> list[pl.Expr]:
    lo = pl.when(pl.col(a) <= pl.col(b)).then(pl.col(a)).otherwise(pl.col(b))
    hi = pl.when(pl.col(a) <= pl.col(b)).then(pl.col(b)).otherwise(pl.col(a))
    return [lo.alias("_a"), hi.alias("_b")]


def run_template(args: argparse.Namespace) -> None:
    local = args.workdir / Path(TEMPLATE_KEY).name
    timings: dict[str, float] = {}
    log("template: download")
    download_verified(TEMPLATE_KEY, local, args.threads, args.chunk_mib, timings)
    hm3 = pl.read_parquet(args.hapmap3).select("CHR", "POS", *sorted_allele_pair("EA", "NEA"))
    log("template: full parse + HapMap3 join")
    t = time.time()
    parts: list[pl.DataFrame] = []
    total = 0
    with rapidgzip.RapidgzipFile(str(local), parallelization=8) as f:
        for block in line_blocks(f):
            rows = parse_lines(block.text).with_row_index("row_number", offset=block.first_row)
            keyed = rows.with_columns(
                chrom_code_expr("Chrom").alias("CHR"), pl.col("Pos").cast(pl.Int32).alias("POS"),
                *sorted_allele_pair("effectAllele", "otherAllele"),
            )
            parts.append(keyed.join(hm3, on=["CHR", "POS", "_a", "_b"], how="semi").select(INDEX_COLS))
            total += rows.height
    index = pl.concat(parts).with_columns(pl.col("row_number").cast(pl.UInt64))
    timings["parse+join"] = time.time() - t
    out = args.workdir / "shared_variant_index.parquet"
    pq.write_table(
        index.to_arrow().replace_schema_metadata({"total_rows": str(total), "template": TEMPLATE_KEY}),
        out, compression="zstd",
    )
    local.unlink()
    log(f"template: {index.height:,} HapMap3 rows of {total:,}; wrote {out.name} "
        f"({out.stat().st_size / 1e6:.1f} MB)")
    log("timings (s): " + ", ".join(f"{k}={v:.1f}" for k, v in timings.items()))
    log(f"peak RSS {peak_rss_mb():.0f} MB")


# ---------------------------------------------------------------- protein mode: one pass


@dataclass(frozen=True)
class PassResult:
    slim: pl.DataFrame
    posmap: pl.DataFrame
    rows: int
    uncompressed_size: int
    index_bytes: bytes


def one_pass(local: Path, index: pl.DataFrame, total_rows: int, parallelization: int,
             chunk_size: int) -> PassResult:
    targets = index["row_number"].to_numpy().astype(np.int64)
    expected_names = index["Name"]
    posmap: list[tuple[int, int, int, int]] = []
    parts: list[pl.DataFrame] = []
    next_mark, rows, end_offset = 0, 0, 0
    with rapidgzip.RapidgzipFile(str(local), parallelization=parallelization,
                                 chunk_size=chunk_size) as f:
        for block in line_blocks(f):
            new_rows, next_mark = posmap_rows(block, max(next_mark, block.offset))
            posmap.extend(new_rows)
            lo = np.searchsorted(targets, block.first_row)
            hi = np.searchsorted(targets, block.first_row + len(block.line_ends))
            parts.append(parse_lines(gather_lines(block, targets[lo:hi])))
            rows = block.first_row + len(block.line_ends)
            end_offset = block.offset + len(block.text)
        buf = io.BytesIO()
        f.export_index(buf, rapidgzip.IndexFormat.INDEXED_GZIP)
    assert rows == total_rows, f"{rows} rows, template had {total_rows}: row order differs"
    gathered = pl.concat(parts)
    assert gathered.height == index.height
    mismatches = int((gathered["Name"] != expected_names).sum())
    assert mismatches == 0, f"{mismatches} gathered rows do not match the shared index"
    return PassResult(
        slim=gathered.select(pl.col("Beta").cast(pl.Float32), pl.col("SE").cast(pl.Float32)),
        posmap=pl.DataFrame(posmap, orient="row", schema={
            "chrom": pl.Int8, "pos": pl.Int64, "line_offset": pl.UInt64, "row_number": pl.UInt64}),
        rows=rows,
        uncompressed_size=end_offset,
        index_bytes=buf.getvalue(),
    )


def write_slim(slim: pl.DataFrame, out: Path) -> None:
    pq.write_table(slim.to_arrow(), out, compression="zstd", use_byte_stream_split=["Beta", "SE"],
                   use_dictionary=False)


# ---------------------------------------------------------------- protein mode: windows


@dataclass(frozen=True)
class Window:
    chrom: int
    start: int
    end: int


@dataclass(frozen=True)
class UncompressedSpan:
    begin: int
    end: int


def uncompressed_span(posmap: pl.DataFrame, window: Window, eof: int) -> UncompressedSpan:
    """Span of uncompressed text containing every row of window. Posmap rows are in file order
    (chromosomes in text order, positions ascending within one). begin = last row on
    window.chrom with pos <= start, else the row before window.chrom's first row; end = first
    later row past window.chrom's section or with pos > end, else EOF."""
    chroms, poss, offs = (posmap[c].to_list() for c in ("chrom", "pos", "line_offset"))
    on_chrom = [i for i, c in enumerate(chroms) if c == window.chrom]
    assert on_chrom, f"no position-map rows on chromosome {window.chrom}"
    before = [i for i in on_chrom if poss[i] <= window.start]
    b = before[-1] if before else max(on_chrom[0] - 1, 0)
    end = next((offs[i] for i in range(b + 1, len(chroms))
                if chroms[i] != window.chrom or poss[i] > window.end), eof)
    return UncompressedSpan(begin=offs[b], end=end)


@dataclass(frozen=True)
class SeekPoint:
    uncmp: int
    cmp: int


def compressed_span(points: list[SeekPoint], span: UncompressedSpan, size: int) -> tuple[int, int]:
    """Compressed [lo, hi) zran reads to decode span: from the seek point at or before span.begin
    (minus one byte, for a mid-byte point) to the first seek point at or after span.end plus
    one indexed_gzip read buffer (1 MiB)."""
    uncmps = [p.uncmp for p in points]
    i = max(bisect_right(uncmps, span.begin) - 1, 0)
    j = bisect_right(uncmps, span.end)
    lo = max(points[i].cmp - 1, 0)
    hi = size if j >= len(points) else min(points[j].cmp + 1 * MIB, size)
    return lo, hi


class RemoteGzipSource(io.RawIOBase):
    """Seekable read-only view of the remote object for indexed_gzip: reads inside the
    prefetched span come from memory, others cost one range request each (counted)."""

    def __init__(self, obj: RemoteObject, stats: FetchStats) -> None:
        self.obj, self.stats = obj, stats
        self.pos, self.cache_lo, self.cache, self.miss_bytes = 0, 0, b"", 0

    def prefetch(self, lo: int, hi: int, threads: int, piece: int) -> None:
        with ThreadPoolExecutor(threads) as pool:
            parts = list(pool.map(lambda b: fetch_range(self.obj.url, b[0], b[1], self.stats),
                                  chunk_bounds(lo, hi, piece)))
        self.cache_lo, self.cache = lo, b"".join(parts)

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.pos

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        self.pos = {io.SEEK_SET: 0, io.SEEK_CUR: self.pos, io.SEEK_END: self.obj.size}[whence] + offset
        return self.pos

    def read(self, n: int = -1) -> bytes:
        n = self.obj.size - self.pos if n < 0 else min(n, self.obj.size - self.pos)
        if n <= 0:
            return b""
        lo, hi = self.cache_lo, self.cache_lo + len(self.cache)
        if lo <= self.pos and self.pos + n <= hi:
            data = self.cache[self.pos - lo : self.pos - lo + n]
        else:
            data = fetch_range(self.obj.url, self.pos, self.pos + n - 1, self.stats)
            self.miss_bytes += n
        self.pos += len(data)
        return data

    def readinto(self, b) -> int:  # noqa: ANN001
        data = self.read(len(b))
        b[: len(data)] = data
        return len(data)


@dataclass(frozen=True)
class WindowReadResult:
    rows: pl.DataFrame
    seconds: float
    stats: FetchStats
    compressed_span: tuple[int, int]
    miss_bytes: int


def remote_window_read(obj: RemoteObject, index_file: Path, posmap: pl.DataFrame,
                       points: list[SeekPoint], eof: int, window: Window, threads: int) -> WindowReadResult:
    t0 = time.time()
    stats = FetchStats()
    span = uncompressed_span(posmap, window, eof)
    lo, hi = compressed_span(points, span, obj.size)
    source = RemoteGzipSource(obj, stats)
    source.prefetch(lo, hi, threads=threads, piece=512 * 1024)
    f = indexed_gzip.IndexedGzipFile(fileobj=source, index_file=str(index_file),
                                     auto_build=False, readbuf_size=1 * MIB)
    f.seek(span.begin)
    text = f.read(span.end - span.begin)
    f.close()
    rows = parse_lines(text).filter(chrom_code_expr("Chrom") == window.chrom,
                                    pl.col("Pos").is_between(window.start, window.end))
    return WindowReadResult(rows=rows, seconds=time.time() - t0, stats=stats,
                            compressed_span=(lo, hi), miss_bytes=source.miss_bytes)


def brute_force_window(local: Path, window: Window) -> pl.DataFrame:
    chrom = "chrX" if window.chrom == 23 else f"chr{window.chrom}"
    awk = f'$1=="{chrom}" && $2>={window.start} && $2<={window.end}'
    out = subprocess.run(f"gzip -dc {local} | awk -F'\\t' '{awk}'", shell=True, check=True,
                         capture_output=True).stdout
    return parse_lines(out)


# ---------------------------------------------------------------- protein mode


def run_protein(args: argparse.Namespace) -> None:
    name = Path(PROTEIN_KEY).name.removesuffix(".txt.gz")
    local = args.workdir / Path(PROTEIN_KEY).name
    index_path = args.workdir / "shared_variant_index.parquet"
    index_tbl = pq.read_table(index_path)
    total_rows = int(index_tbl.schema.metadata[b"total_rows"])
    index = pl.from_arrow(index_tbl)
    timings: dict[str, float] = {}
    t_build = time.time()

    log("protein: stage 1-2 presign, download, md5")
    obj = download_verified(PROTEIN_KEY, local, args.threads, args.chunk_mib, timings)

    log(f"protein: stage 3 one rapidgzip pass (parallelization {args.parallelization}, "
        f"chunk_size {args.chunk_size_mib} MiB)")
    t = time.time()
    res = one_pass(local, index, total_rows, args.parallelization, args.chunk_size_mib * MIB)
    timings["one pass"] = time.time() - t
    log(f"  {res.rows:,} rows ({res.uncompressed_size / 1e6:.0f} MB) checked against the shared "
        f"index; gathered {res.slim.height:,} rows; posmap {res.posmap.height} rows; "
        f"peak RSS so far {peak_rss_mb():.0f} MB")

    log("protein: stage 4 write outputs")
    t = time.time()
    slim_path = args.workdir / f"{name}.slim.parquet"
    write_slim(res.slim, slim_path)
    index_raw = args.workdir / f"{name}.gzidx"
    index_raw.write_bytes(res.index_bytes)
    index_zst = args.workdir / f"{name}.gzidx.zst"
    index_zst.write_bytes(zstandard.ZstdCompressor(level=3).compress(res.index_bytes))
    posmap_path = args.workdir / f"{name}.posmap.parquet"
    res.posmap.write_parquet(posmap_path, compression="zstd")
    timings["write"] = time.time() - t
    build_seconds = time.time() - t_build
    for p in (slim_path, index_raw, index_zst, posmap_path):
        log(f"  {p.name}: {p.stat().st_size / 1e6:.2f} MB")
    log(f"BUILD per protein: {build_seconds:.1f}s; "
        + ", ".join(f"{k}={v:.1f}" for k, v in timings.items())
        + f"; peak RSS {peak_rss_mb():.0f} MB")

    with indexed_gzip.IndexedGzipFile(str(local), index_file=str(index_raw)) as f:
        points = [SeekPoint(uncmp=u, cmp=c) for u, c in f.raw.seek_points()]
    log(f"  seek points: {len(points)}")

    z = (pl.Series(res.slim["Beta"]) / pl.Series(res.slim["SE"])).abs()
    lead = index.row(int(z.arg_max()), named=True)
    windows = {
        "lead": Window(chrom_code(lead["Chrom"]), max(lead["Pos"] - HALF_WINDOW_BP, 1),
                       lead["Pos"] + HALF_WINDOW_BP),
        "late-file chr9": Window(9, 100_000_000, 101_000_000),
    }
    log(f"  lead index variant {lead['Name']} |z|={float(z.max()):.1f}")
    ok = True
    for label, window in windows.items():
        log(f"protein: stage 5 remote window read [{label}] {window}")
        w = remote_window_read(obj, index_raw, res.posmap, points, res.uncompressed_size,
                               window, threads=args.window_threads)
        lo, hi = w.compressed_span
        log(f"  {w.rows.height:,} rows in {w.seconds:.1f}s; fetched {w.stats.bytes / 1e6:.1f} MB "
            f"in {w.stats.requests} requests (outside prefetch {w.miss_bytes / 1e6:.2f} MB); "
            f"compressed span {(hi - lo) / 1e6:.1f} MB")
        expected = brute_force_window(local, window)
        same = expected.equals(w.rows)
        ok &= same
        log(f"protein: stage 6 brute force [{label}]: {expected.height:,} rows; identical: {same}")
    local.unlink()
    log(f"peak RSS (whole process) {peak_rss_mb():.0f} MB")
    log("ALL WINDOWS MATCH" if ok else "WINDOW MISMATCH")
    sys.exit(0 if ok else 1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["template", "protein"])
    ap.add_argument("--workdir", type=Path, required=True)
    ap.add_argument("--hapmap3", type=Path)
    ap.add_argument("--threads", type=int, default=128)
    ap.add_argument("--chunk-mib", type=int, default=4)
    ap.add_argument("--parallelization", type=int, default=8)
    ap.add_argument("--chunk-size-mib", type=int, default=16)
    ap.add_argument("--window-threads", type=int, default=16)
    args = ap.parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    if args.mode == "template":
        assert args.hapmap3 is not None, "template mode needs --hapmap3"
        run_template(args)
    else:
        run_protein(args)


if __name__ == "__main__":
    main()
