"""
Download one file from a presigned URL with many parallel HTTP Range requests, using only the
standard library (threads + urllib). Each worker writes its chunk at its offset in a
preallocated output file. Reports throughput and the md5 of the result.

Usage: python -I python_parallel_range_download.py <presigned_url> <out_path> <n_threads>
"""

import hashlib
import os
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

CHUNK = 8 * 1024 * 1024


def total_size(url: str) -> int:
    req = urllib.request.Request(url, headers={"Range": "bytes=0-0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers["Content-Range"].split("/")[1])


def fetch(url: str, fd: int, start: int, end: int) -> int:
    for attempt in range(8):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(req, timeout=120) as r:
                assert r.status == 206, r.status
                data = r.read()
            assert len(data) == end - start + 1, (len(data), start, end)
            os.pwrite(fd, data, start)
            return len(data)
        except (OSError, AssertionError) as e:
            print(f"retry {attempt} for {start}: {e}", file=sys.stderr)
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"chunk at {start} failed")


def main() -> None:
    url, out, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
    size = total_size(url)
    fd = os.open(out, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    os.ftruncate(fd, size)
    t0 = time.time()
    with ThreadPoolExecutor(n) as pool:
        done = sum(
            pool.map(
                lambda s: fetch(url, fd, s, min(s + CHUNK, size) - 1),
                range(0, size, CHUNK),
            )
        )
    dt = time.time() - t0
    os.close(fd)
    md5 = hashlib.md5()
    with open(out, "rb") as f:
        while b := f.read(1 << 24):
            md5.update(b)
    print(f"python threads={n}: {done / 1e6:.0f} MB in {dt:.0f}s = {done / 1e6 / dt:.1f} MB/s md5={md5.hexdigest()}")


if __name__ == "__main__":
    main()
