"""
Fetch dense protein summary statistics for a few locus windows, then discard the download.

The DenseWindowFetcher protocol is what the colocalization Task depends on, so any proteomics
source (UKB-PPP, CSF pQTL, ...) can be plugged in. Fetched frames are in the common window
schema (WINDOW_SCHEMA, positions in GRCh38) but NOT yet oriented to the reference genome:
EA/NEA are whatever the source reports, and orientation happens downstream.

PppWindowFetcher downloads a protein's UKB-PPP tar from Synapse once per call, reads only the
chromosome files the windows need (one at a time, only the needed columns), filters each
window, and deletes the tar and the extracted files before returning. In PPP's regenie files
ALLELE1 is the effect allele and A1FREQ its frequency; ALLELE0 is the other allele.
"""

import shutil
from collections.abc import Sequence
from pathlib import Path
from typing import Protocol

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.ppp_database.build_slim_protein_parquet_task import (
    PppProteinFile,
    extract_protein_tar,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_CHROM_NAME_FOR_CODE,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.protein_coloc_constants import WINDOW_SCHEMA
from mecfs_bio.constants.regenie_constants import (
    REGENIE_A1FREQ_COL,
    REGENIE_ALLELE0_COL,
    REGENIE_ALLELE1_COL,
    REGENIE_BETA_COL,
    REGENIE_CHROM_COL,
    REGENIE_GENPOS_COL,
    REGENIE_N_COL,
    REGENIE_SE_COL,
)

_MAX_CHROM = 23
_REGENIE_TO_WINDOW = {
    REGENIE_CHROM_COL: GWASLAB_CHROM_COL,
    REGENIE_GENPOS_COL: GWASLAB_POS_COL,
    REGENIE_ALLELE1_COL: GWASLAB_EFFECT_ALLELE_COL,
    REGENIE_ALLELE0_COL: GWASLAB_NON_EFFECT_ALLELE_COL,
    REGENIE_BETA_COL: GWASLAB_BETA_COL,
    REGENIE_SE_COL: GWASLAB_SE_COL,
    REGENIE_A1FREQ_COL: GWASLAB_EFFECT_ALLELE_FREQ_COL,
    REGENIE_N_COL: GWASLAB_SAMPLE_SIZE_COLUMN,
}


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


class DenseWindowFetcher(Protocol):
    def fetch(
        self,
        protein_id: str,
        windows: Sequence[LocusWindow],
        scratch_dir: Path,
        wf: WF,
    ) -> dict[str, pl.DataFrame]:
        """Return one WINDOW_SCHEMA frame per window, keyed by locus_id."""
        ...


def _chromosome_file(extracted_dir: Path, chrom: int) -> Path:
    name = GWASLAB_CHROM_NAME_FOR_CODE.get(chrom, None) or str(chrom)
    matches = list(extracted_dir.rglob(f"*chr{name}_*.gz"))
    assert len(matches) == 1, (
        f"expected exactly one regenie file for chr{name} under {extracted_dir}, "
        f"found {matches}"
    )
    return matches[0]


def _read_chromosome(path: Path) -> pl.DataFrame:
    return (
        pl.read_csv(path, separator=" ", columns=list(_REGENIE_TO_WINDOW))
        .rename(_REGENIE_TO_WINDOW)
        .cast(WINDOW_SCHEMA)
        .select(WINDOW_SCHEMA.names())
    )


def _slice_window(chromosome_rows: pl.DataFrame, window: LocusWindow) -> pl.DataFrame:
    return chromosome_rows.filter(
        pl.col(GWASLAB_CHROM_COL) == window.chrom,
        pl.col(GWASLAB_POS_COL).is_between(window.start, window.end),
    )


@frozen(slots=True)
class PppWindowFetcher:
    """
    DenseWindowFetcher for UKB-PPP.

    proteins: the manifest rows of every protein that may be requested, looked up by OID.
    """

    proteins: tuple[PppProteinFile, ...]

    def _protein(self, protein_id: str) -> PppProteinFile:
        matches = [p for p in self.proteins if p.oid == protein_id]
        assert len(matches) == 1, f"expected one manifest row for {protein_id}"
        return matches[0]

    def fetch(
        self,
        protein_id: str,
        windows: Sequence[LocusWindow],
        scratch_dir: Path,
        wf: WF,
    ) -> dict[str, pl.DataFrame]:
        protein = self._protein(protein_id)
        download_dir = scratch_dir / f"dense_{protein_id}"
        download_dir.mkdir(parents=True, exist_ok=True)
        try:
            tar_path = wf.download_from_synapse(
                protein.synid, download_dir, protein.tar_filename
            )
            extracted_dir = extract_protein_tar(tar_path, download_dir)
            tar_path.unlink()
            frames: dict[str, pl.DataFrame] = {}
            for chrom in sorted({w.chrom for w in windows}):
                rows = _read_chromosome(_chromosome_file(extracted_dir, chrom))
                for window in windows:
                    if window.chrom == chrom:
                        frames[window.locus_id] = _slice_window(rows, window)
            return frames
        finally:
            shutil.rmtree(download_dir, ignore_errors=True)
