"""
Source-specific readers turning a gene-coordinate table into the map protein id -> GRCh38
gene interval that candidate selection uses to call a pair cis or trans.

A protein with no entry has no cis window; its pairs get cis status "unknown".
"""

from typing import Protocol

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.ppp_ldsc.gene_coords import read_gene_coords
from mecfs_bio.constants.protein_coloc_constants import ProteinId
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval


class GeneCoordsReader(Protocol):
    def __call__(self, table: pl.DataFrame) -> dict[ProteinId, GenomicInterval]: ...


@frozen(slots=True)
class St3GeneCoordsReader:
    """UKB-PPP: the extracted Sun et al. 2023 ST3 sheet, keyed by Olink id."""

    def __call__(self, table: pl.DataFrame) -> dict[ProteinId, GenomicInterval]:
        return {
            ProteinId(str(oid)): interval
            for oid, interval in read_gene_coords(table).items()
        }
