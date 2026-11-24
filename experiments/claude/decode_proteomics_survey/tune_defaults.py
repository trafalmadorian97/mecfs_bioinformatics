"""
Tuning measurements for the deCODE defaults (plan Task 17, step 2):

  pixi r python experiments/claude/decode_proteomics_survey/tune_defaults.py chunk <scratch_dir>
      download AP4M1 with 128 connections at chunk_bytes 2, 4 and 8 MiB (twice each);
  pixi r python experiments/claude/decode_proteomics_survey/tune_defaults.py build
      build the first 24 SMP proteins of the manifest (needed for the parallelism test);
  pixi r python experiments/claude/decode_proteomics_survey/tune_defaults.py parallel
      read one 1 Mb window from each of the 24 through ParallelWindowFetcher at
      max_parallel 4, 6 and 8 (twice each).
"""

import sys
import tempfile
import time
from contextlib import closing
from pathlib import Path

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_database_protein_files import (
    HAPMAP_3_DECODE_SMP_DATABASE,
)
from mecfs_bio.assets.gwas.decode_somascan.raw.decode_variant_info import (
    DECODE_VARIANT_INFO,
)
from mecfs_bio.build_system.task.protein_coloc.decode_window_reader import (
    DecodeProteinEntry,
    DecodeWindowReader,
)
from mecfs_bio.build_system.task.protein_coloc.dense_window_fetch import (
    ParallelWindowFetcher,
    ProteinWindowRequest,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

MIB = 1024 * 1024
AP4M1 = "10076-1"
N_PROTEINS = 24
WINDOW = GenomicInterval(11, 116_300_000, 117_300_000)


def _log(msg: str) -> None:
    print(msg, flush=True)


def chunk(scratch: Path) -> None:
    (task,) = [
        t for t in HAPMAP_3_DECODE_SMP_DATABASE.protein_tasks if t.protein.seq_id == AP4M1
    ]
    wf = make_wf()
    protein = task.protein
    for _ in range(2):
        for chunk_mib in (2, 4, 8):
            out = scratch / f"ap4m1_{chunk_mib}.txt.gz"
            t0 = time.perf_counter()
            wf.download_url_parallel(
                url=wf.s3_access.presigned_url(protein.ref, 3600),
                size=protein.size,
                local_path=out,
                md5_hash=protein.md5,
                n_connections=128,
                chunk_bytes=chunk_mib * MIB,
            )
            seconds = time.perf_counter() - t0
            _log(
                f"chunk {chunk_mib} MiB: {seconds:.1f}s "
                f"({protein.size / MIB / seconds:.0f} MiB/s)"
            )
            out.unlink()


def _tasks():  # noqa: ANN202
    return HAPMAP_3_DECODE_SMP_DATABASE.protein_tasks[:N_PROTEINS]


def build() -> None:
    t0 = time.perf_counter()
    DEFAULT_RUNNER.run(list(_tasks()), incremental_save=True)
    _log(f"built {N_PROTEINS} proteins in {time.perf_counter() - t0:.0f}s")


def parallel() -> None:
    tasks = _tasks()
    assets = DEFAULT_RUNNER.run([*tasks, DECODE_VARIANT_INFO])
    reader = DecodeWindowReader(
        proteins=tuple(DecodeProteinEntry(t) for t in tasks),
        variant_info_task=DECODE_VARIANT_INFO,
    )
    requests = [
        ProteinWindowRequest(protein_id=str(t.protein.seq_id), windows=(WINDOW,))
        for t in tasks
    ]
    wf = make_wf()
    for _ in range(2):
        for max_parallel in (4, 6, 8):
            fetcher = ParallelWindowFetcher(reader=reader, max_parallel=max_parallel)
            with tempfile.TemporaryDirectory() as scratch:
                t0 = time.perf_counter()
                with closing(
                    fetcher.fetch_all(requests, Path(scratch), lambda a: assets[a], wf)
                ) as results:
                    rows = sum(f.frames[WINDOW].height for f in results)
                _log(
                    f"max_parallel {max_parallel}: {time.perf_counter() - t0:.1f}s for "
                    f"{N_PROTEINS} proteins ({rows} rows)"
                )


def main() -> None:
    mode = sys.argv[1]
    if mode == "chunk":
        chunk(Path(sys.argv[2]))
    elif mode == "build":
        build()
    else:
        parallel()


if __name__ == "__main__":
    main()
