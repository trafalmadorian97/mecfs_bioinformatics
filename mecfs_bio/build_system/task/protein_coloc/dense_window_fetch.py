"""
The interface between the colocalization Task and a proteomics GWAS source.

The Task hands a DenseWindowFetcher every selected protein with its locus windows, and reads
back one result per protein, in request order. Fetched frames are in the common window schema
(WINDOW_SCHEMA, positions in GRCh38) but NOT yet oriented to the reference genome: EA/NEA are
whatever the source reports, and orientation happens downstream.

How a source gets its data (concurrency, retries) is its own business. Sources that fetch one
protein at a time through a ProteinWindowSource can use fetch_proteins_in_parallel, which runs
up to max_parallel of those single-protein fetches at once, ahead of the consumer, and still
yields results in request order. A fetch mostly waits on the network, so running several at
once overlaps the waits.
"""

import time
from collections.abc import Generator, Sequence
from concurrent.futures import Future, ThreadPoolExecutor
from functools import partial
from pathlib import Path
from typing import Protocol

import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.util.retry.call_with_retries import call_with_retries

logger = structlog.get_logger()

_MAX_CHROM = 23
# Fetches hit remote services; network and filesystem errors are worth retrying.
FETCH_TRANSIENT_ERRORS: tuple[type[Exception], ...] = (OSError,)


@frozen(slots=True)
class LocusWindow:
    """A closed GRCh38 interval [start, end] on chromosome chrom (X coded as 23)."""

    locus_id: str
    chrom: int
    start: int
    end: int

    def __attrs_post_init__(self) -> None:
        assert self.start <= self.end, f"window start > end: {self}"
        assert 1 <= self.chrom <= _MAX_CHROM, f"chromosome out of range: {self}"


@frozen(slots=True)
class ProteinWindowRequest:
    """One protein and the locus windows to fetch for it."""

    protein_id: str
    windows: tuple[LocusWindow, ...]


@frozen(slots=True)
class ProteinWindowFrames:
    """The fetched windows of one request: one WINDOW_SCHEMA frame per window, keyed by
    locus_id."""

    request: ProteinWindowRequest
    frames: dict[str, pl.DataFrame]


class DenseWindowFetcher(Protocol):
    def fetch_all(
        self,
        requests: Sequence[ProteinWindowRequest],
        scratch_dir: Path,
        wf: WF,
    ) -> Generator[ProteinWindowFrames]:
        """
        One result per request, in request order.

        Close the generator when stopping early (e.g. with contextlib.closing), so that work
        on requests not yet consumed stops.
        """
        ...


class ProteinWindowSource(Protocol):
    def fetch_protein(
        self,
        protein_id: str,
        windows: Sequence[LocusWindow],
        scratch_dir: Path,
        wf: WF,
    ) -> dict[str, pl.DataFrame]:
        """Return one WINDOW_SCHEMA frame per window, keyed by locus_id."""
        ...


class SleepFn(Protocol):
    def __call__(self, seconds: float, /) -> None: ...


def _fetch_with_retries(
    source: ProteinWindowSource,
    request: ProteinWindowRequest,
    scratch_dir: Path,
    wf: WF,
    sleep: SleepFn,
) -> ProteinWindowFrames:
    frames = call_with_retries(
        partial(
            source.fetch_protein, request.protein_id, request.windows, scratch_dir, wf
        ),
        retry_on=FETCH_TRANSIENT_ERRORS,
        sleep=sleep,
    )
    logger.info(
        "fetched dense protein windows",
        protein_id=request.protein_id,
        n_windows=len(request.windows),
    )
    return ProteinWindowFrames(request=request, frames=frames)


def fetch_proteins_in_parallel(
    source: ProteinWindowSource,
    requests: Sequence[ProteinWindowRequest],
    scratch_dir: Path,
    wf: WF,
    max_parallel: int,
    sleep: SleepFn = time.sleep,
) -> Generator[ProteinWindowFrames]:
    """
    Fetch every request through source, up to max_parallel at once, retrying transient
    errors (backing off with sleep), and yield the results in request order.

    All fetches are queued at once and run ahead of the consumer. Fetched windows are small,
    so letting fetches outpace a slow consumer costs little memory. A fetch that still fails
    after its retries raises when its result is reached. Closing the generator cancels the
    fetches not yet started and waits for those running.
    """
    assert max_parallel >= 1, f"max_parallel must be >= 1, got {max_parallel}"
    protein_ids = [r.protein_id for r in requests]
    assert len(set(protein_ids)) == len(protein_ids), "duplicate protein in requests"
    pool = ThreadPoolExecutor(max_workers=max_parallel)
    try:
        futures: list[Future[ProteinWindowFrames]] = [
            pool.submit(_fetch_with_retries, source, request, scratch_dir, wf, sleep)
            for request in requests
        ]
        for future in futures:
            yield future.result()
    finally:
        pool.shutdown(wait=True, cancel_futures=True)
