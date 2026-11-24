"""
Gene coordinates for the Western et al. 2024 CSF aptamers: each aptamer's NCBI Entrez gene id
(from the SomaScan aptamer table) and gene symbol, looked up in the GRCh38 MAGMA gene-location
table by the rules of magma_gene_intervals.

An aptamer with no interval has cis status "unknown" in the coloc results.
"""

import re
from typing import Final

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.csf_database.build_slim_aptamer_parquet_task import (
    CsfAptamerFile,
)
from mecfs_bio.build_system.task.protein_coloc.candidate_selection import (
    SelectionConfig,
)
from mecfs_bio.build_system.task.protein_coloc.magma_gene_intervals import (
    GeneKeys,
    gene_interval,
    magma_genes,
)
from mecfs_bio.constants.protein_coloc_constants import ProteinId
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

# Multi-gene ids in the aptamer table are separated by "|" or a space.
_ID_SEPARATOR: Final = re.compile(r"[| ]+")


def _gene_keys(aptamer: CsfAptamerFile) -> GeneKeys:
    ids = (
        frozenset()
        if aptamer.entrez_gene_id is None
        else frozenset(_ID_SEPARATOR.split(aptamer.entrez_gene_id.strip()))
    )
    return GeneKeys(entrez_ids=ids, symbol=aptamer.entrez_gene_symbol)


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
        genes = magma_genes(table)
        intervals = {
            ProteinId(aptamer.seq_id): gene_interval(
                _gene_keys(aptamer), genes, self.conflict_distance_bp
            )
            for aptamer in self.aptamers
        }
        return {pid: iv for pid, iv in intervals.items() if iv is not None}
