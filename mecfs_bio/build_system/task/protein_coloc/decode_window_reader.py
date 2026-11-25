"""
ProteinWindowReader for the deCODE SomaScan proteomics GWAS: reads trait-locus windows of one
protein's remote sumstats file by range requests, using the seek index and position map its
BuildDecodeSomascanProteinTask stored. Wrap it in a ParallelWindowFetcher,
DECODE_MAX_PARALLEL_FETCHES at once, to fetch many proteins.

For each window, in file order:
  1. the position map gives an uncompressed span holding every row of the window. Chromosomes
     are in text order (chr1, chr10, ..., chr9, chrX), so the span is found by comparing
     (text-order rank, position) keys, never by chromosome code alone;
  2. the seek index gives the compressed byte range that decodes that span. The range is
     prefetched with parallel range requests, then indexed_gzip decodes the span from the
     prefetched bytes (any read outside them costs one more range request);
  3. the lines are parsed and filtered to the window. Rows whose other allele is the
     missing-allele marker are dropped, and so are rows that are not clean (keep_clean_rows):
     at multi-allelic and overlapping-record positions deCODE tests one allele against the
     pooled rest, which a trait GWAS row does not measure, so such rows must not enter coloc.

EAF is the true effect allele frequency from the clean variant_info entries (resolve_effect_
allele_frequency), null where the rule does not determine it; deCODE's ImpMAF is never used as
EAF.

Before reading, the object's size and md5 are checked against the position map's metadata. A
changed object raises RemoteObjectChanged, which is not an OSError, so ParallelWindowFetcher
fails at once rather than retrying: the stored index no longer describes the file.
"""

import io
from abc import ABC, abstractmethod
from bisect import bisect_right
from collections.abc import Buffer, Sequence
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import indexed_gzip
import polars as pl
import pyarrow.parquet as pq
from attrs import frozen

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.decode_database.build_decode_somascan_protein_task import (
    DECODE_POSMAP_FILENAME,
    DECODE_SEEK_INDEX_FILENAME,
    POSMAP_COMPRESSED_SIZE_KEY,
    POSMAP_MD5_KEY,
    BuildDecodeSomascanProteinTask,
)
from mecfs_bio.build_system.task.decode_database.build_decode_variant_info_task import (
    keep_clean_rows,
)
from mecfs_bio.build_system.task.decode_database.decode_sumstats_format import (
    chrom_code_expr,
    parse_decode_lines,
)
from mecfs_bio.build_system.task.decode_database.effect_allele_frequency import (
    resolve_effect_allele_frequency,
)
from mecfs_bio.build_system.task.decode_database.line_blocks import MIB
from mecfs_bio.build_system.task.decode_database.one_pass import (
    POSMAP_CHROM_COL,
    POSMAP_LINE_OFFSET_COL,
    POSMAP_POS_COL,
)
from mecfs_bio.build_system.task.protein_coloc.dense_window_fetch import (
    ProteinWindowReader,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.csf_database_constants import SeqId
from mecfs_bio.constants.decode_somascan_constants import (
    DECODE_BETA_COL,
    DECODE_CHR_X_CODE,
    DECODE_CHROM_COL,
    DECODE_EFFECT_ALLELE_COL,
    DECODE_IMP_MAF_COL,
    DECODE_MISSING_ALLELE,
    DECODE_N_COL,
    DECODE_OTHER_ALLELE_COL,
    DECODE_POS_COL,
    DECODE_SE_COL,
    VI_CHR_COL,
    VI_POS_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.protein_coloc_constants import WINDOW_SCHEMA
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval
from mecfs_bio.util.download.parallel_ranges import (
    REQUESTS_RANGE_GETTER,
    RangeGetter,
    chunk_bounds,
    fetch_range_with_retries,
)
from mecfs_bio.util.file_io.gzip_seek_index import SeekPoint, load_seek_index

# Proteins fetched at once. One 1 Mb window from each of 24 proteins: 4 at once 19-21 s, 6 13-17 s,
# 8 12-13 s, 12 9-12 s, 16 11-13 s; peak RSS 1.7 GB.
DECODE_MAX_PARALLEL_FETCHES = 8
# indexed_gzip's compressed read buffer: decoding a span may read this far past its end point.
_READBUF_BYTES = 1 * MIB
_PRESIGN_EXPIRES_S = 3600

# Rank of each chromosome code in the files' text order (chr1, chr10, ..., chr9, chrX).
_TEXT_ORDER_RANK: dict[int, int] = {
    code: rank
    for rank, code in enumerate(
        sorted(
            range(1, DECODE_CHR_X_CODE + 1),
            key=lambda c: "chrX" if c == DECODE_CHR_X_CODE else f"chr{c}",
        )
    )
}


class RemoteObjectChanged(Exception):
    """The remote file no longer matches the one the stored seek index was built from."""


class RangeSource(ABC):
    @abstractmethod
    def read_range(self, start: int, length: int) -> bytes:
        """length bytes of the object from offset start."""
        ...


class RangeSourceFactory(ABC):
    @abstractmethod
    def open(self, url: str, size: int) -> RangeSource:
        """A RangeSource over the object of size bytes at url."""
        ...


@frozen(slots=True)
class _PresignedRangeSource(RangeSource):
    url: str
    getter: RangeGetter

    def read_range(self, start: int, length: int) -> bytes:
        return fetch_range_with_retries(
            self.getter, self.url, start, start + length - 1
        )


@frozen(slots=True)
class PresignedRangeSourceFactory(RangeSourceFactory):
    """Range GETs of a presigned URL, retrying transient failures."""

    getter: RangeGetter = REQUESTS_RANGE_GETTER

    def open(self, url: str, size: int) -> RangeSource:
        return _PresignedRangeSource(url=url, getter=self.getter)


@frozen(slots=True)
class UncompressedSpan:
    """
    begin: uncompressed offset of a line start at or before the window's first row.
    end: uncompressed offset at or after the end of the window's last row; None = end of file.
    """

    begin: int
    end: int | None


@frozen(slots=True)
class CompressedSpan:
    """Compressed bytes [lo, hi) that indexed_gzip reads to decode an UncompressedSpan."""

    lo: int
    hi: int


def _window_key(chrom: int, pos: int) -> tuple[int, int]:
    return _TEXT_ORDER_RANK[chrom], pos


def uncompressed_span(
    posmap: pl.DataFrame, window: GenomicInterval
) -> UncompressedSpan:
    """
    posmap rows are line starts in file order, the first at the first data line. begin is the
    last row whose (rank, pos) key is below (rank, window.start): every row at or after
    window.start lies after it, even when several rows share a position. end is the first row
    whose key is above (rank, window.end), else the end of the file.
    """
    keys = [
        _window_key(c, p)
        for c, p in zip(
            posmap[POSMAP_CHROM_COL].to_list(),
            posmap[POSMAP_POS_COL].to_list(),
            strict=True,
        )
    ]
    offsets = posmap[POSMAP_LINE_OFFSET_COL].to_list()
    start_key = _window_key(window.chrom, window.start)
    end_key = _window_key(window.chrom, window.end)
    before = [i for i, k in enumerate(keys) if k < start_key]
    b = before[-1] if before else 0
    end = next((offsets[i] for i in range(b + 1, len(keys)) if keys[i] > end_key), None)
    return UncompressedSpan(begin=offsets[b], end=end)


def compressed_span(
    points: Sequence[SeekPoint], span: UncompressedSpan, size: int
) -> CompressedSpan:
    """
    From the seek point at or before span.begin (one byte early, as a point may start
    mid-byte) to the first seek point past span.end plus one indexed_gzip read buffer, or to
    the end of the object of size bytes.
    """
    uncompressed = [p.uncompressed for p in points]
    i = max(bisect_right(uncompressed, span.begin) - 1, 0)
    lo = max(points[i].compressed - 1, 0)
    if span.end is None:
        return CompressedSpan(lo=lo, hi=size)
    j = bisect_right(uncompressed, span.end)
    hi = size if j >= len(points) else min(points[j].compressed + _READBUF_BYTES, size)
    return CompressedSpan(lo=lo, hi=hi)


class RemoteGzipSource(io.RawIOBase):
    """
    Seekable read-only view of a remote object for indexed_gzip. Reads inside the prefetched
    span come from memory; any other read costs one range request.
    """

    def __init__(self, source: RangeSource, size: int) -> None:
        super().__init__()
        self._source = source
        self._size = size
        self._pos = 0
        self._cache_lo = 0
        self._cache = b""

    def prefetch(
        self, span: CompressedSpan, pool: ThreadPoolExecutor, piece_bytes: int
    ) -> None:
        """
        Fetch [span.lo, span.hi) in pieces of piece_bytes on pool's threads. Reusing one pool
        across windows reuses each thread's HTTP connection.
        """
        pieces = chunk_bounds(span.lo, span.hi, piece_bytes)
        parts = list(
            pool.map(lambda c: self._source.read_range(c.start, c.length), pieces)
        )
        self._cache_lo, self._cache = span.lo, b"".join(parts)

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self._pos

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        base = {io.SEEK_SET: 0, io.SEEK_CUR: self._pos, io.SEEK_END: self._size}[whence]
        self._pos = base + offset
        return self._pos

    def read(self, size: int | None = -1) -> bytes:
        remaining = self._size - self._pos
        n = remaining if size is None or size < 0 else min(size, remaining)
        if n <= 0:
            return b""
        cache_hi = self._cache_lo + len(self._cache)
        if self._cache_lo <= self._pos and self._pos + n <= cache_hi:
            start = self._pos - self._cache_lo
            data = self._cache[start : start + n]
        else:
            data = self._source.read_range(self._pos, n)
        self._pos += len(data)
        return data

    def readinto(self, buffer: Buffer, /) -> int:
        view = memoryview(buffer).cast("B")
        data = self.read(len(view))
        view[: len(data)] = data
        return len(data)


@frozen(slots=True)
class DecodeProteinEntry:
    """One protein the reader may be asked for: the task that stored its seek index."""

    task: BuildDecodeSomascanProteinTask

    @property
    def seq_id(self) -> SeqId:
        return self.task.protein.seq_id


@frozen(slots=True)
class _StoredIndex:
    posmap: pl.DataFrame
    seek_index: bytes
    size: int
    md5: str


def _load_stored_index(asset: DirectoryAsset) -> _StoredIndex:
    posmap_path = asset.path / DECODE_POSMAP_FILENAME
    metadata = pq.read_schema(posmap_path).metadata
    assert metadata is not None, f"{posmap_path} has no metadata"
    return _StoredIndex(
        posmap=pl.read_parquet(posmap_path),
        seek_index=(asset.path / DECODE_SEEK_INDEX_FILENAME).read_bytes(),
        size=int(metadata[POSMAP_COMPRESSED_SIZE_KEY.encode()]),
        md5=metadata[POSMAP_MD5_KEY.encode()].decode(),
    )


def _check_unchanged(
    task: BuildDecodeSomascanProteinTask, stored: _StoredIndex, wf: WF
) -> None:
    info = wf.s3_access.head(task.protein.ref)
    if info.size != stored.size or info.md5 != stored.md5:
        raise RemoteObjectChanged(
            f"{task.protein.key} changed since its seek index was built: size {info.size} "
            f"md5 {info.md5}, indexed size {stored.size} md5 {stored.md5}"
        )


def _read_text(f: indexed_gzip.IndexedGzipFile, span: UncompressedSpan) -> bytes:
    f.seek(span.begin)
    return f.read() if span.end is None else f.read(span.end - span.begin)


def _window_rows(text: bytes, window: GenomicInterval) -> pl.DataFrame:
    return parse_decode_lines(text).filter(
        (chrom_code_expr(DECODE_CHROM_COL) == window.chrom)
        & pl.col(DECODE_POS_COL).is_between(window.start, window.end)
        & (pl.col(DECODE_OTHER_ALLELE_COL) != DECODE_MISSING_ALLELE)
    )


def _variant_info_rows(
    variant_info: pl.LazyFrame, window: GenomicInterval
) -> pl.DataFrame:
    return variant_info.filter(
        (pl.col(VI_CHR_COL) == window.chrom)
        & pl.col(VI_POS_COL).is_between(window.start, window.end)
    ).collect()


def _to_window_frame(rows: pl.DataFrame, variant_info: pl.DataFrame) -> pl.DataFrame:
    frame = rows.select(
        chrom_code_expr(DECODE_CHROM_COL).alias(GWASLAB_CHROM_COL),
        pl.col(DECODE_POS_COL).alias(GWASLAB_POS_COL),
        pl.col(DECODE_EFFECT_ALLELE_COL).alias(GWASLAB_EFFECT_ALLELE_COL),
        pl.col(DECODE_OTHER_ALLELE_COL).alias(GWASLAB_NON_EFFECT_ALLELE_COL),
        pl.col(DECODE_BETA_COL).alias(GWASLAB_BETA_COL),
        pl.col(DECODE_SE_COL).alias(GWASLAB_SE_COL),
        pl.col(DECODE_N_COL).alias(GWASLAB_SAMPLE_SIZE_COLUMN),
        pl.col(DECODE_IMP_MAF_COL),
    )
    frame = keep_clean_rows(frame, variant_info)
    eaf = resolve_effect_allele_frequency(frame, variant_info)
    return (
        frame.with_columns(eaf.alias(GWASLAB_EFFECT_ALLELE_FREQ_COL))
        .cast(WINDOW_SCHEMA)
        .select(WINDOW_SCHEMA.names())
    )


def _in_file_order(windows: Sequence[GenomicInterval]) -> list[GenomicInterval]:
    return sorted(windows, key=lambda w: (_window_key(w.chrom, w.start), w.end))


@frozen(slots=True)
class DecodeWindowReader(ProteinWindowReader):
    """
    proteins: every protein that may be requested, looked up by seq id (the coloc protein id).
    clean_variant_info_task: the BuildDecodeCleanVariantInfoTask whose entries decide which
        rows are clean and give their true EAFs.
    range_source_factory: how range requests reach the presigned URL.
    prefetch_connections, prefetch_piece_bytes: parallel range requests per window prefetch,
        and the bytes each requests.
    """

    proteins: tuple[DecodeProteinEntry, ...]
    clean_variant_info_task: Task
    range_source_factory: RangeSourceFactory = PresignedRangeSourceFactory()
    prefetch_connections: int = 16
    prefetch_piece_bytes: int = 512 * 1024

    @property
    def deps(self) -> list[Task]:
        return [entry.task for entry in self.proteins] + [self.clean_variant_info_task]

    def _task(self, protein_id: str) -> BuildDecodeSomascanProteinTask:
        matches = [e.task for e in self.proteins if e.seq_id == protein_id]
        assert len(matches) == 1, f"expected one deCODE protein for {protein_id}"
        return matches[0]

    def fetch_protein(
        self,
        protein_id: str,
        windows: Sequence[GenomicInterval],
        scratch_dir: Path,
        fetch: Fetch,
        wf: WF,
    ) -> dict[GenomicInterval, pl.DataFrame]:
        """One protein's windows, keyed by window."""
        task = self._task(protein_id)
        asset = fetch(task.asset_id)
        assert isinstance(asset, DirectoryAsset)
        stored = _load_stored_index(asset)
        _check_unchanged(task, stored, wf)
        source = RemoteGzipSource(
            self.range_source_factory.open(
                wf.s3_access.presigned_url(task.protein.ref, _PRESIGN_EXPIRES_S),
                stored.size,
            ),
            stored.size,
        )
        variant_info = scan_dataframe_asset(
            fetch(self.clean_variant_info_task.asset_id),
            meta=self.clean_variant_info_task.meta,
            parquet_backend="polars",
        ).to_native()
        assert isinstance(variant_info, pl.LazyFrame)
        f = indexed_gzip.IndexedGzipFile(
            fileobj=source, auto_build=False, readbuf_size=_READBUF_BYTES
        )
        try:
            f.import_index(fileobj=load_seek_index(stored.seek_index))
            points = [
                SeekPoint(uncompressed=u, compressed=c) for u, c in f.raw.seek_points()
            ]
            with ThreadPoolExecutor(self.prefetch_connections) as pool:
                return {
                    window: self._read_window(
                        f, source, pool, stored, points, window, variant_info
                    )
                    for window in _in_file_order(windows)
                }
        finally:
            f.close()

    def _read_window(
        self,
        f: indexed_gzip.IndexedGzipFile,
        source: RemoteGzipSource,
        pool: ThreadPoolExecutor,
        stored: _StoredIndex,
        points: list[SeekPoint],
        window: GenomicInterval,
        variant_info: pl.LazyFrame,
    ) -> pl.DataFrame:
        span = uncompressed_span(stored.posmap, window)
        source.prefetch(
            compressed_span(points, span, stored.size),
            pool=pool,
            piece_bytes=self.prefetch_piece_bytes,
        )
        rows = _window_rows(_read_text(f, span), window)
        return _to_window_frame(rows, _variant_info_rows(variant_info, window))
