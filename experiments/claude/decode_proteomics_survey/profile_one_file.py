"""
Profile one deCODE proteomics sumstats file streamed on stdin (plain gzip, sorted by Chrom/Pos).

Reports rows, uncompressed size, throughput, and the compressed byte offset (to within one
read chunk) at which each chromosome starts: the cost of reaching a window by decompressing
from the start of the file.

Usage: aws --endpoint-url https://s3-ext.decode.is s3 cp s3://largescaleplasma-2023/<key> - \
    | python profile_one_file.py
"""

import sys
import time
import zlib

CHUNK = 1 << 20


def main() -> None:
    d = zlib.decompressobj(31)
    t0 = time.time()
    comp = 0
    uncomp = 0
    rows = 0
    tail = b""
    current = None
    starts: list[tuple[str, int, int]] = []
    ns: set[bytes] = set()
    bang = 0
    while chunk := sys.stdin.buffer.read(CHUNK):
        comp += len(chunk)
        out = tail + d.decompress(chunk)
        uncomp += len(out) - len(tail)
        lines = out.split(b"\n")
        tail = lines.pop()
        for line in lines:
            if line.startswith(b"Chrom"):
                continue
            f = line.split(b"\t")
            rows += 1
            if f[0] != current:
                current = f[0]
                starts.append((current.decode(), comp, rows))
            if rows % 1000 == 0:
                ns.add(f[10])
            if f[5] == b"!":
                bang += 1
    dt = time.time() - t0
    print(f"rows={rows} compressed={comp / 1e6:.0f}MB uncompressed={uncomp / 1e6:.0f}MB "
          f"seconds={dt:.0f} MB/s(compressed)={comp / 1e6 / dt:.1f} multiallelic_bang={bang}")
    print(f"sampled distinct N: {sorted(ns)[:10]} (count {len(ns)})")
    for chrom, off, r in starts:
        print(f"{chrom}\tcompressed_offset_MB={off / 1e6:.0f}\tfrac={off / comp:.3f}\trow={r}")


if __name__ == "__main__":
    main()
