"""
The GRCh38 interval of a measured protein's gene, from its NCBI Entrez gene id(s) and gene
symbol, looked up in the MAGMA gene-location table (keyed by Entrez id). Shared by every
proteomics source that labels its proteins with Entrez ids.

Rules, per protein:
1. The id is in MAGMA: use that gene's interval, unless the protein's gene symbol names a
   different MAGMA gene on another chromosome or more than conflict_distance_bp away. Then the
   two keys disagree about where the target is (for example CBS, whose id points at a
   duplicated copy 36 Mb from the canonical locus), and the protein gets no interval.
   Readthroughs and adjacent paralogs, close together, keep the id's interval.
2. The id is not in MAGMA (retired or replaced ids): fall back to the symbol, if it names
   exactly one MAGMA gene.
3. Several distinct ids (complexes, ambiguous paralog sets), no id, or neither key found:
   no interval.
"""

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.ppp_ldsc.gene_coords import parse_chrom
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_ENTREZ_ID_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
)
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval


@frozen(slots=True)
class GeneKeys:
    """entrez_ids: the protein's distinct Entrez gene ids (may be empty); symbol: its gene
    symbol."""

    entrez_ids: frozenset[str]
    symbol: str


@frozen(slots=True)
class MagmaGenes:
    """MAGMA genes by Entrez id, and by symbol for symbols naming exactly one gene."""

    by_id: dict[str, GenomicInterval]
    by_symbol: dict[str, GenomicInterval]


def magma_genes(table: pl.DataFrame) -> MagmaGenes:
    """Index the MAGMA gene-location table; genes on unplaced contigs are skipped."""
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
    return MagmaGenes(by_id=by_id, by_symbol=by_symbol)


def _far_apart(a: GenomicInterval, b: GenomicInterval, distance_bp: int) -> bool:
    if a.chrom != b.chrom:
        return True
    gap = max(a.start, b.start) - min(a.end, b.end)
    return gap > distance_bp


def gene_interval(
    keys: GeneKeys, genes: MagmaGenes, conflict_distance_bp: int
) -> GenomicInterval | None:
    """The protein's gene interval by the rules above, or None."""
    if len(keys.entrez_ids) > 1:
        return None
    by_id = genes.by_id.get(next(iter(keys.entrez_ids))) if keys.entrez_ids else None
    by_symbol = genes.by_symbol.get(keys.symbol)
    if by_id is None:
        return by_symbol
    if by_symbol is not None and _far_apart(by_id, by_symbol, conflict_distance_bp):
        return None
    return by_id
