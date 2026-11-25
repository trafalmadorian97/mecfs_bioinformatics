"""
Checks of the smoke-run outputs against brute-force parses of a local copy of the AP4M1 file
(plan Task 17, steps 1, 3 and 4):

  pixi r python experiments/claude/decode_proteomics_survey/verify_smoke.py <AP4M1.txt.gz>

1. slim == Beta/SE/N of the file's row with each index NAME (assert_frame_equal);
   significant_rows row count and size.
2. three DecodeWindowReader windows == an awk scan of the file (EAF aside).
3. the spread of N (p05, p95) over the slim.
"""

import io
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import polars as pl
from polars.testing import assert_frame_equal

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_database_protein_files import (
    HAPMAP_3_DECODE_SMP_DATABASE,
)
from mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_index import (
    HAPMAP_3_DECODE_SMP_INDEX,
)
from mecfs_bio.assets.gwas.decode_somascan.raw.decode_variant_info import (
    DECODE_VARIANT_INFO,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.decode_database.build_decode_somascan_protein_task import (
    DECODE_SIGNIFICANT_ROWS_FILENAME,
    DECODE_SLIM_FILENAME,
)
from mecfs_bio.build_system.task.decode_database.decode_sumstats_format import (
    parse_decode_lines,
)
from mecfs_bio.build_system.task.protein_coloc.decode_window_reader import (
    DecodeProteinEntry,
    DecodeWindowReader,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.decode_somascan_constants import DECODE_INDEX_NAME_COL
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
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

SEQ_ID = "10076-1"
WINDOWS = (
    GenomicInterval(1, 1, 500_000),
    GenomicInterval(11, 116_500_000, 117_000_000),
    GenomicInterval(23, 154_000_000, 156_100_000),
)


def _log(msg: str) -> None:
    print(msg, flush=True)


def check_slim(local: Path, index_path: Path, protein_dir: Path) -> None:
    t0 = time.perf_counter()
    cols = subprocess.run(
        f"gzip -dc {local} | cut -f3,7,10,11",
        shell=True,
        check=True,
        capture_output=True,
    ).stdout
    full = pl.read_csv(io.BytesIO(cols), separator="\t", null_values=["."])
    index = pl.read_parquet(index_path, columns=[DECODE_INDEX_NAME_COL])
    expected = index.join(
        full,
        left_on=DECODE_INDEX_NAME_COL,
        right_on="Name",
        how="left",
        maintain_order="left",
    ).select(
        pl.col("Beta").cast(pl.Float32).alias(GWASLAB_BETA_COL),
        pl.col("SE").cast(pl.Float32).alias(GWASLAB_SE_COL),
        pl.col("N").cast(pl.Float32).alias(GWASLAB_SAMPLE_SIZE_COLUMN),
    )
    slim = pl.read_parquet(protein_dir / DECODE_SLIM_FILENAME)
    assert_frame_equal(slim, expected)
    n = slim[GWASLAB_SAMPLE_SIZE_COLUMN].to_numpy().astype(float)
    _log(
        f"slim == brute force ({slim.height} rows, {time.perf_counter() - t0:.0f}s); "
        f"file rows {full.height}; N p05 {np.quantile(n, 0.05):.0f} median "
        f"{np.median(n):.0f} p95 {np.quantile(n, 0.95):.0f} min {n.min():.0f} max {n.max():.0f}"
    )
    sig = protein_dir / DECODE_SIGNIFICANT_ROWS_FILENAME
    _log(
        f"significant_rows: {pl.read_parquet(sig).height} rows, {sig.stat().st_size} bytes"
    )


def _awk_window(local: Path, window: GenomicInterval) -> pl.DataFrame:
    chrom = "chrX" if window.chrom == 23 else f"chr{window.chrom}"
    awk = f'$1=="{chrom}" && $2>={window.start} && $2<={window.end} && $6!="!"'
    out = subprocess.run(
        f"gzip -dc {local} | awk -F'\\t' '{awk}'",
        shell=True,
        check=True,
        capture_output=True,
    ).stdout
    rows = parse_decode_lines(out)
    return rows.select(
        pl.lit(window.chrom, dtype=pl.Int64).alias(GWASLAB_CHROM_COL),
        pl.col("Pos").alias(GWASLAB_POS_COL),
        pl.col("effectAllele").alias(GWASLAB_EFFECT_ALLELE_COL),
        pl.col("otherAllele").alias(GWASLAB_NON_EFFECT_ALLELE_COL),
        pl.col("Beta").alias(GWASLAB_BETA_COL),
        pl.col("SE").alias(GWASLAB_SE_COL),
        pl.col("N").alias(GWASLAB_SAMPLE_SIZE_COLUMN),
    )


def check_windows(local: Path, protein_task, assets) -> None:  # noqa: ANN001
    reader = DecodeWindowReader(
        proteins=(DecodeProteinEntry(protein_task),),
        variant_info_task=DECODE_VARIANT_INFO,
    )
    with tempfile.TemporaryDirectory() as scratch:
        t0 = time.perf_counter()
        frames = reader.fetch_protein(
            SEQ_ID, WINDOWS, Path(scratch), lambda asset_id: assets[asset_id], make_wf()
        )
        _log(f"read {len(WINDOWS)} remote windows in {time.perf_counter() - t0:.1f}s")
    for window in WINDOWS:
        got = frames[window]
        assert_frame_equal(
            got.drop(GWASLAB_EFFECT_ALLELE_FREQ_COL), _awk_window(local, window)
        )
        eaf = got[GWASLAB_EFFECT_ALLELE_FREQ_COL]
        _log(f"window {window}: {got.height} rows == awk; EAF null {eaf.null_count()}")


def main() -> None:
    local = Path(sys.argv[1])
    (protein_task,) = [
        t
        for t in HAPMAP_3_DECODE_SMP_DATABASE.protein_tasks
        if t.protein.seq_id == SEQ_ID
    ]
    assets = DEFAULT_RUNNER.run(
        [protein_task, HAPMAP_3_DECODE_SMP_INDEX, DECODE_VARIANT_INFO]
    )
    index = assets[HAPMAP_3_DECODE_SMP_INDEX.asset_id]
    protein = assets[protein_task.asset_id]
    assert isinstance(index, FileAsset) and isinstance(protein, DirectoryAsset)
    _log(f"index rows: {pl.read_parquet(index.path, columns=[GWASLAB_POS_COL]).height}")
    check_slim(local, index.path, protein.path)
    check_windows(local, protein_task, assets)


if __name__ == "__main__":
    main()
