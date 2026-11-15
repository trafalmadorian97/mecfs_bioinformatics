"""
Check live that PppWindowReader's ranged reads give the same windows as the full download.

Uses AMY1A_AMY1B_AMY1C, whose tar stores member names in GNU long-name entries. Fetches two
windows on two chromosomes with the ranged fetcher, then downloads and extracts the whole tar
and slices the same windows from the extracted chromosome files, and compares. Then fetches
the same windows for AMY1A_AMY1B_AMY1C and five more proteins through fetch_all (4 at once)
and checks the results come back in request order and equal fetch_protein's, one at a time.

Run: pixi r python experiments/claude/protein_coloc/verify_ppp_ranged_window_fetch.py 2>&1 | tee experiments/claude/protein_coloc/verify_ppp_ranged_window_fetch.log
"""

import tempfile
import time
from pathlib import Path

import polars as pl
from polars.testing import assert_frame_equal

from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval
from mecfs_bio.asset_generator.ukbb_ppp_slim_protein_asset_generator import (
    EUR_DISCOVERY_PPP_MANIFEST_PATH,
    load_ppp_protein_files,
)
from mecfs_bio.build_system.task.ppp_database.build_slim_protein_parquet_task import (
    PppProteinFile,
    extract_protein_tar,
)
from mecfs_bio.build_system.task.protein_coloc.dense_window_fetch import (
    ParallelWindowFetcher,
    ProteinWindowRequest,
)
from mecfs_bio.build_system.task.protein_coloc.ppp_window_reader import (
    _REGENIE_TO_WINDOW,
    PPP_MAX_PARALLEL_FETCHES,
    PppWindowReader,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import GWASLAB_POS_COL
from mecfs_bio.constants.protein_coloc_constants import WINDOW_SCHEMA

_GENE_PREFIX = "AMY1A_AMY1B_AMY1C"
_WINDOWS = (
    GenomicInterval(chrom=1, start=103_500_000, end=104_500_000),
    GenomicInterval(chrom=22, start=30_000_000, end=31_000_000),
)


def _full_download_windows(
    protein_synid: str, tar_name: str, scratch: Path
) -> dict[GenomicInterval, pl.DataFrame]:
    tar_path = make_wf().download_from_synapse(protein_synid, scratch, tar_name)
    extracted = extract_protein_tar(tar_path, scratch)
    frames = {}
    for window in _WINDOWS:
        (path,) = list(extracted.rglob(f"*chr{window.chrom}_*.gz"))
        rows = (
            pl.read_csv(path, separator=" ", columns=list(_REGENIE_TO_WINDOW))
            .rename(_REGENIE_TO_WINDOW)
            .cast(WINDOW_SCHEMA)
            .select(WINDOW_SCHEMA.names())
        )
        frames[window] = rows.filter(
            pl.col(GWASLAB_POS_COL).is_between(window.start, window.end)
        )
    return frames


def main() -> None:
    proteins = load_ppp_protein_files(EUR_DISCOVERY_PPP_MANIFEST_PATH)
    (protein,) = [p for p in proteins if p.tar_filename.startswith(_GENE_PREFIX)]
    reader = PppWindowReader(proteins=(protein,))
    with tempfile.TemporaryDirectory() as tmp:
        scratch = Path(tmp)
        start = time.monotonic()
        ranged = reader.fetch_protein(protein.oid, _WINDOWS, scratch, make_wf())
        print(f"ranged fetch: {time.monotonic() - start:.1f}s")
        start = time.monotonic()
        full = _full_download_windows(protein.synid, protein.tar_filename, scratch)
        print(f"full download + extract: {time.monotonic() - start:.1f}s")
    for window in _WINDOWS:
        print(f"{window}: {ranged[window].height} rows")
        assert ranged[window].height > 0
        assert_frame_equal(ranged[window], full[window])
    print("ranged and full-download windows are identical")
    _check_fetch_all(proteins, protein)


_N_EXTRA_PROTEINS = 5


def _check_fetch_all(
    proteins: tuple[PppProteinFile, ...], amy1: PppProteinFile
) -> None:
    batch = (amy1, *proteins[20 : 20 + _N_EXTRA_PROTEINS])
    reader = PppWindowReader(proteins=batch)
    fetcher = ParallelWindowFetcher(reader=reader, max_parallel=PPP_MAX_PARALLEL_FETCHES)
    requests = [
        ProteinWindowRequest(protein_id=p.oid, windows=_WINDOWS) for p in batch
    ]
    with tempfile.TemporaryDirectory() as tmp:
        scratch = Path(tmp)
        start = time.monotonic()
        parallel = list(fetcher.fetch_all(requests, scratch, make_wf()))
        print(f"fetch_all, {len(batch)} proteins, 4 at once: {time.monotonic() - start:.1f}s")
        start = time.monotonic()
        serial = [
            reader.fetch_protein(r.protein_id, r.windows, scratch, make_wf())
            for r in requests
        ]
        print(f"fetch_protein one at a time: {time.monotonic() - start:.1f}s")
    assert [r.request for r in parallel] == requests
    for result, frames in zip(parallel, serial, strict=True):
        for window in _WINDOWS:
            assert_frame_equal(result.frames[window], frames[window])
    print("fetch_all results are in request order and equal fetch_protein's")


if __name__ == "__main__":
    main()
