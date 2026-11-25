"""
ProteinWindowReader for UKB-PPP: fetch one protein's dense summary statistics for a few locus
windows. Wrap it in a ParallelWindowFetcher, PPP_MAX_PARALLEL_FETCHES at once, to fetch many
proteins.

PppWindowReader never downloads a protein's whole UKB-PPP tar (about 550 MB, one gzipped
regenie file per chromosome). It lists the tar's members with one small HTTP Range read per
member header, then range-reads only the chromosome files the windows need, one at a time,
decompressing each in memory and keeping only the needed columns. A chromosome file is 1-8%
of the tar. The files are plain gzip rather than BGZF, so a chromosome file is always read
whole. In PPP's regenie files ALLELE1 is the effect allele and A1FREQ its frequency; ALLELE0
is the other allele.
"""

import fnmatch
import gzip
import io
from collections.abc import Sequence
from pathlib import Path, PurePosixPath

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.ppp_database.build_slim_protein_parquet_task import (
    PppProteinFile,
)
from mecfs_bio.build_system.task.protein_coloc.dense_window_fetch import (
    ProteinWindowReader,
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
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval
from mecfs_bio.util.file_io.ranged_tar import (
    RangeReader,
    TarMemberExtent,
    list_tar_members,
    read_tar_member,
)

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


def _chromosome_member(
    members: Sequence[TarMemberExtent], chrom: int
) -> TarMemberExtent:
    name = GWASLAB_CHROM_NAME_FOR_CODE.get(chrom, None) or str(chrom)
    pattern = f"*chr{name}_*.gz"
    matches = [
        m for m in members if fnmatch.fnmatch(PurePosixPath(m.name).name, pattern)
    ]
    assert len(matches) == 1, (
        f"expected exactly one regenie file for chr{name} in the tar, found {matches}"
    )
    return matches[0]


def _read_chromosome(reader: RangeReader, member: TarMemberExtent) -> pl.DataFrame:
    text = gzip.decompress(read_tar_member(reader, member))
    return (
        pl.read_csv(io.BytesIO(text), separator=" ", columns=list(_REGENIE_TO_WINDOW))
        .rename(_REGENIE_TO_WINDOW)
        .cast(WINDOW_SCHEMA)
        .select(WINDOW_SCHEMA.names())
    )


def _slice_window(
    chromosome_rows: pl.DataFrame, window: GenomicInterval
) -> pl.DataFrame:
    return chromosome_rows.filter(
        pl.col(GWASLAB_CHROM_COL) == window.chrom,
        pl.col(GWASLAB_POS_COL).is_between(window.start, window.end),
    )


def _slice_windows_on_chromosome(
    chromosome_rows: pl.DataFrame, windows: Sequence[GenomicInterval], chrom: int
) -> dict[GenomicInterval, pl.DataFrame]:
    return {w: _slice_window(chromosome_rows, w) for w in windows if w.chrom == chrom}


# How many PPP proteins to fetch at once. Timed on 12-protein batches, 2, 4 and 8 at once were
# 2.0, 3.4 and 4.6 times faster than one at a time, with 8 near CPU-bound (gzip decompression
# and CSV parsing) and peak memory about 1 GB.
PPP_MAX_PARALLEL_FETCHES = 4


@frozen(slots=True)
class PppWindowReader(ProteinWindowReader):
    """
    ProteinWindowReader for UKB-PPP.

    proteins: the manifest rows of every protein that may be requested, looked up by OID.
    Nothing is written to scratch_dir.
    """

    proteins: tuple[PppProteinFile, ...]

    def _protein(self, protein_id: str) -> PppProteinFile:
        matches = [p for p in self.proteins if p.oid == protein_id]
        assert len(matches) == 1, f"expected one manifest row for {protein_id}"
        return matches[0]

    def fetch_protein(
        self,
        protein_id: str,
        windows: Sequence[GenomicInterval],
        scratch_dir: Path,
        wf: WF,
    ) -> dict[GenomicInterval, pl.DataFrame]:
        """One protein's windows, keyed by window."""
        protein = self._protein(protein_id)
        reader = wf.open_synapse_range_reader(protein.synid)
        members = list_tar_members(reader)
        frames: dict[GenomicInterval, pl.DataFrame] = {}
        for chrom in sorted({w.chrom for w in windows}):
            rows = _read_chromosome(reader, _chromosome_member(members, chrom))
            frames.update(_slice_windows_on_chromosome(rows, windows, chrom))
        return frames
