"""
Check live that PppWindowFetcher's ranged reads give the same windows as the full download.

Uses AMY1A_AMY1B_AMY1C, whose tar stores member names in GNU long-name entries. Fetches two
windows on two chromosomes with the ranged fetcher, then downloads and extracts the whole tar
and slices the same windows from the extracted chromosome files, and compares.

Run: pixi r python experiments/claude/protein_coloc/verify_ppp_ranged_window_fetch.py 2>&1 | tee experiments/claude/protein_coloc/verify_ppp_ranged_window_fetch.log
"""

import tempfile
import time
from pathlib import Path

import polars as pl
from polars.testing import assert_frame_equal

from mecfs_bio.asset_generator.ukbb_ppp_slim_protein_asset_generator import (
    EUR_DISCOVERY_PPP_MANIFEST_PATH,
    load_ppp_protein_files,
)
from mecfs_bio.build_system.task.ppp_database.build_slim_protein_parquet_task import (
    extract_protein_tar,
)
from mecfs_bio.build_system.task.protein_coloc.ppp_window_fetcher import (
    _REGENIE_TO_WINDOW,
    LocusWindow,
    PppWindowFetcher,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import GWASLAB_POS_COL
from mecfs_bio.constants.protein_coloc_constants import WINDOW_SCHEMA

_GENE_PREFIX = "AMY1A_AMY1B_AMY1C"
_WINDOWS = (
    LocusWindow(locus_id="chr1_amy", chrom=1, start=103_500_000, end=104_500_000),
    LocusWindow(locus_id="chr22_any", chrom=22, start=30_000_000, end=31_000_000),
)


def _full_download_windows(
    protein_synid: str, tar_name: str, scratch: Path
) -> dict[str, pl.DataFrame]:
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
        frames[window.locus_id] = rows.filter(
            pl.col(GWASLAB_POS_COL).is_between(window.start, window.end)
        )
    return frames


def main() -> None:
    proteins = load_ppp_protein_files(EUR_DISCOVERY_PPP_MANIFEST_PATH)
    (protein,) = [p for p in proteins if p.tar_filename.startswith(_GENE_PREFIX)]
    fetcher = PppWindowFetcher(proteins=(protein,))
    with tempfile.TemporaryDirectory() as tmp:
        scratch = Path(tmp)
        start = time.monotonic()
        ranged = fetcher.fetch(protein.oid, _WINDOWS, scratch, make_wf())
        print(f"ranged fetch: {time.monotonic() - start:.1f}s")
        start = time.monotonic()
        full = _full_download_windows(protein.synid, protein.tar_filename, scratch)
        print(f"full download + extract: {time.monotonic() - start:.1f}s")
    for window in _WINDOWS:
        print(f"{window.locus_id}: {ranged[window.locus_id].height} rows")
        assert ranged[window.locus_id].height > 0
        assert_frame_equal(ranged[window.locus_id], full[window.locus_id])
    print("ranged and full-download windows are identical")


if __name__ == "__main__":
    main()
