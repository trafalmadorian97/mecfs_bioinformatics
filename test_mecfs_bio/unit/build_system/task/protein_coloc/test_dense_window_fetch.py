import threading
from collections.abc import Sequence
from pathlib import Path

import polars as pl
import pytest

from mecfs_bio.build_system.task.protein_coloc.dense_window_fetch import (
    ProteinWindowRequest,
    fetch_proteins_in_parallel,
)
from mecfs_bio.build_system.wf.base_wf import WF, make_wf
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

_BARRIER_TIMEOUT_SECONDS = 10
_GATE_RELEASE_SECONDS = 0.2
_WINDOW = GenomicInterval(chrom=1, start=1, end=100)


def _requests(n: int) -> list[ProteinWindowRequest]:
    return [
        ProteinWindowRequest(protein_id=f"P{i}", windows=(_WINDOW,)) for i in range(n)
    ]


def _frames_for(protein_id: str) -> dict[GenomicInterval, pl.DataFrame]:
    return {_WINDOW: pl.DataFrame({"protein": [protein_id]})}


class _FakeSource:
    """
    barrier: when given, every fetch waits on it, so fetches complete only if that many run
        at once.
    release_order: when given, fetches finish in this protein order, whatever order they
        started in.
    n_failures_before_success: each protein's first fetches raise OSError this many times.
    """

    def __init__(
        self,
        barrier: threading.Barrier | None = None,
        release_order: Sequence[str] | None = None,
        n_failures_before_success: int = 0,
    ) -> None:
        self.barrier = barrier
        self.release_order = list(release_order) if release_order else None
        self.n_failures_before_success = n_failures_before_success
        self.calls: list[str] = []
        self._condition = threading.Condition()

    def fetch_protein(
        self,
        protein_id: str,
        windows: Sequence[GenomicInterval],
        scratch_dir: Path,
        wf: WF,
    ) -> dict[GenomicInterval, pl.DataFrame]:
        with self._condition:
            self.calls.append(protein_id)
            n_calls = self.calls.count(protein_id)
        if n_calls <= self.n_failures_before_success:
            raise OSError("transient")
        if self.barrier is not None:
            self.barrier.wait()
        release_order = self.release_order
        if release_order is not None:
            with self._condition:
                assert self._condition.wait_for(
                    lambda: release_order[0] == protein_id,
                    timeout=_BARRIER_TIMEOUT_SECONDS,
                )
                release_order.pop(0)
                self._condition.notify_all()
        return _frames_for(protein_id)


def _no_sleep(seconds: float) -> None:
    pass


def _fetch_all(
    source: _FakeSource, requests: list[ProteinWindowRequest], max_parallel: int
) -> list[str]:
    results = fetch_proteins_in_parallel(
        source,
        requests,
        Path("unused"),
        make_wf(),
        max_parallel=max_parallel,
        sleep=_no_sleep,
    )
    return [r.request.protein_id for r in results]


def test_fetches_run_at_once_up_to_max_parallel() -> None:
    requests = _requests(3)
    source = _FakeSource(barrier=threading.Barrier(3, timeout=_BARRIER_TIMEOUT_SECONDS))
    assert _fetch_all(source, requests, max_parallel=3) == ["P0", "P1", "P2"]


def test_results_come_in_request_order_when_fetches_finish_out_of_order() -> None:
    requests = _requests(3)
    source = _FakeSource(release_order=["P2", "P0", "P1"])
    results = list(
        fetch_proteins_in_parallel(
            source, requests, Path("unused"), make_wf(), max_parallel=3
        )
    )
    assert [r.request for r in results] == requests
    assert [r.frames[_WINDOW]["protein"].item() for r in results] == ["P0", "P1", "P2"]


def test_transient_errors_are_retried() -> None:
    source = _FakeSource(n_failures_before_success=1)
    assert _fetch_all(source, _requests(2), max_parallel=1) == ["P0", "P1"]
    assert sorted(source.calls) == ["P0", "P0", "P1", "P1"]


class _GatedSource:
    """Fetches of gated_protein block until gate is set; other fetches return at once."""

    def __init__(self, gated_protein: str) -> None:
        self.gated_protein = gated_protein
        self.gate = threading.Event()
        self.calls: list[str] = []

    def fetch_protein(
        self,
        protein_id: str,
        windows: Sequence[GenomicInterval],
        scratch_dir: Path,
        wf: WF,
    ) -> dict[GenomicInterval, pl.DataFrame]:
        self.calls.append(protein_id)
        if protein_id == self.gated_protein:
            assert self.gate.wait(timeout=_BARRIER_TIMEOUT_SECONDS)
        return _frames_for(protein_id)


def test_closing_early_cancels_fetches_not_yet_started() -> None:
    # One worker: P1 runs (blocked on the gate) while P2-P4 wait in the queue.
    source = _GatedSource(gated_protein="P1")
    results = fetch_proteins_in_parallel(
        source, _requests(5), Path("unused"), make_wf(), max_parallel=1
    )
    assert next(results).request.protein_id == "P0"
    # Release P1 shortly after close() has cancelled the queue, so close() can finish.
    release = threading.Timer(_GATE_RELEASE_SECONDS, source.gate.set)
    release.start()
    results.close()
    release.join()
    assert source.calls == ["P0", "P1"]


def test_duplicate_proteins_are_rejected() -> None:
    requests = _requests(1) * 2
    with pytest.raises(AssertionError):
        _fetch_all(_FakeSource(), requests, max_parallel=1)
