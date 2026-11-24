"""
Seek indexes for plain (single-stream) gzip files, built by rapidgzip and read by indexed_gzip.

rapidgzip decompresses a gzip file in parallel and can export the seek points it found as a
GZIDX v1 index (the indexed_gzip format). indexed_gzip imports that index and can then
decompress from any seek point onward, so a reader can extract a byte span of the
uncompressed text while reading only the compressed bytes around it.

A seek point pairs an uncompressed offset with the compressed offset at which decoding can
resume (plus the 32 KiB window stored in the index). The index is stored zstd-compressed,
since the windows compress well.
"""

import io
from pathlib import Path
from typing import BinaryIO

import indexed_gzip
import rapidgzip
import zstandard
from attrs import frozen

_ZSTD_LEVEL = 3


@frozen(slots=True)
class SeekPoint:
    """
    uncompressed: offset into the decompressed text.
    compressed: offset into the gzip file where decoding of that text can resume.
    """

    uncompressed: int
    compressed: int


def export_index_from_open_file(f: rapidgzip.RapidgzipFile) -> bytes:
    """The zstd-compressed GZIDX index of a RapidgzipFile that has been read to the end."""
    buf = io.BytesIO()
    f.export_index(buf, rapidgzip.IndexFormat.INDEXED_GZIP)
    return zstandard.ZstdCompressor(level=_ZSTD_LEVEL).compress(buf.getvalue())


def export_seek_index(
    gz_path: Path, out_path: Path, parallelization: int, chunk_size: int
) -> None:
    """
    Decompress gz_path once with rapidgzip and write its zstd-compressed GZIDX index.
    chunk_size is the compressed distance in bytes between seek points.
    """
    with rapidgzip.RapidgzipFile(
        str(gz_path), parallelization=parallelization, chunk_size=chunk_size
    ) as f:
        while f.read(16 * 1024 * 1024):
            pass
        out_path.write_bytes(export_index_from_open_file(f))


def load_seek_index(blob: bytes) -> io.BytesIO:
    """A file object over the decompressed GZIDX index, for IndexedGzipFile.import_index."""
    return io.BytesIO(zstandard.ZstdDecompressor().decompress(blob))


def seek_points(blob: bytes, source: BinaryIO) -> tuple[SeekPoint, ...]:
    """The seek points of the zstd-compressed GZIDX index blob, for the gzip data in source."""
    f = indexed_gzip.IndexedGzipFile(fileobj=source, auto_build=False)
    try:
        f.import_index(fileobj=load_seek_index(blob))
        return tuple(
            SeekPoint(uncompressed=u, compressed=c) for u, c in f.raw.seek_points()
        )
    finally:
        f.close()
