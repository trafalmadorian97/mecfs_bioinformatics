"""
Gene coordinates for the Western et al. 2024 CSF aptamers: each aptamer's NCBI Entrez gene id
(from the SomaScan aptamer table) joined to the GRCh38 MAGMA gene-location table, which is
keyed by Entrez id.

Rules, per aptamer:
1. The id is in MAGMA: use that gene's interval, unless the aptamer's gene symbol names a
   different MAGMA gene on another chromosome or more than conflict_distance_bp away. Then the
   two keys disagree about where the target is (for example CBS, whose id points at a
   duplicated copy 36 Mb from the canonical locus), and the aptamer gets no interval.
   Readthroughs and adjacent paralogs, close together, keep the id's interval.
2. The id is not in MAGMA (retired or replaced ids): fall back to the symbol, if it names
   exactly one MAGMA gene.
3. Several distinct ids (complexes, ambiguous paralog sets), no id, or neither key found:
   no interval.

An aptamer with no interval has cis status "unknown" in the coloc results.
"""

import re
from typing import Final

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.csf_database.build_slim_aptamer_parquet_task import (
    CsfAptamerFile,
)
from mecfs_bio.build_system.task.ppp_ldsc.gene_coords import parse_chrom
from mecfs_bio.build_system.task.protein_coloc.candidate_selection import (
    SelectionConfig,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_ENTREZ_ID_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
)
from mecfs_bio.constants.protein_coloc_constants import ProteinId
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

# Multi-gene ids in the aptamer table are separated by "|" or a space.
_ID_SEPARATOR: Final = re.compile(r"[| ]+")


@frozen(slots=True)
class _MagmaGenes:
    by_id: dict[str, GenomicInterval]
    by_symbol: dict[str, GenomicInterval]


def _magma_genes(table: pl.DataFrame) -> _MagmaGenes:
    by_id: dict[str, GenomicInterval] = {}
    symbol_rows: dict[str, list[GenomicInterval]] = {}
    for row in table.select(
        pl.col(GENE_INFO_ENTREZ_ID_COL).cast(pl.String),
        pl.col(GENE_INFO_CHROM_COL).cast(pl.String),
        pl.col(GENE_INFO_START_COL).cast(pl.Int64),
        pl.col(GENE_INFO_END_COL).cast(pl.Int64),
        pl.col(GENE_INFO_NAME_COL).cast(pl.String),
    ).iter_rows(named=True):
        chrom = parse_chrom(row[GENE_INFO_CHROM_COL])
        if chrom is None:
            continue
        interval = GenomicInterval(
            chrom=chrom, start=row[GENE_INFO_START_COL], end=row[GENE_INFO_END_COL]
        )
        by_id[row[GENE_INFO_ENTREZ_ID_COL]] = interval
        symbol_rows.setdefault(row[GENE_INFO_NAME_COL], []).append(interval)
    by_symbol = {s: rows[0] for s, rows in symbol_rows.items() if len(rows) == 1}
    return _MagmaGenes(by_id=by_id, by_symbol=by_symbol)


def _entrez_ids(aptamer: CsfAptamerFile) -> set[str]:
    if aptamer.entrez_gene_id is None:
        return set()
    return set(_ID_SEPARATOR.split(aptamer.entrez_gene_id.strip()))


def _far_apart(a: GenomicInterval, b: GenomicInterval, distance_bp: int) -> bool:
    if a.chrom != b.chrom:
        return True
    gap = max(a.start, b.start) - min(a.end, b.end)
    return gap > distance_bp


def _aptamer_interval(
    aptamer: CsfAptamerFile, genes: _MagmaGenes, conflict_distance_bp: int
) -> GenomicInterval | None:
    ids = _entrez_ids(aptamer)
    if len(ids) > 1:
        return None
    by_id = genes.by_id.get(next(iter(ids))) if ids else None
    by_symbol = genes.by_symbol.get(aptamer.entrez_gene_symbol)
    if by_id is None:
        return by_symbol
    if by_symbol is not None and _far_apart(by_id, by_symbol, conflict_distance_bp):
        return None
    return by_id


@frozen(slots=True)
class CsfGeneCoordsReader:
    """
    aptamers: every aptamer that may appear in the hit table, keyed by seq id.
    conflict_distance_bp: how far the id's and the symbol's genes may lie apart before the
        aptamer is left without an interval; the cis window, so a conflict that could not
        change a cis call is tolerated.
    """

    aptamers: tuple[CsfAptamerFile, ...]
    conflict_distance_bp: int = SelectionConfig().cis_window_bp

    def __call__(self, table: pl.DataFrame) -> dict[ProteinId, GenomicInterval]:
        genes = _magma_genes(table)
        intervals = {
            ProteinId(aptamer.seq_id): _aptamer_interval(
                aptamer, genes, self.conflict_distance_bp
            )
            for aptamer in self.aptamers
        }
        return {pid: iv for pid, iv in intervals.items() if iv is not None}
