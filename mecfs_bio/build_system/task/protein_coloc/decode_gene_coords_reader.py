"""
Gene coordinates for the deCODE SomaScan assays: each assay's NCBI Entrez gene id(s) (from its
UniProt accession via HGNC, in the manifest) and gene symbol (from the file name), looked up in
the GRCh38 MAGMA gene-location table by the rules of magma_gene_intervals.

An assay with no interval has cis status "unknown" in the coloc results.
"""

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.decode_database.decode_somascan_protein import (
    DecodeSomascanProtein,
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

# The manifest joins several Entrez ids with "|".
_ID_SEPARATOR = "|"


def _gene_keys(protein: DecodeSomascanProtein) -> GeneKeys:
    ids = (
        frozenset()
        if protein.entrez_gene_id is None
        else frozenset(protein.entrez_gene_id.split(_ID_SEPARATOR))
    )
    return GeneKeys(entrez_ids=ids, symbol=protein.gene_symbol)


@frozen(slots=True)
class DecodeGeneCoordsReader:
    """
    proteins: every assay that may appear in the hit table, keyed by seq id.
    conflict_distance_bp: how far the id's and the symbol's genes may lie apart before the
        assay is left without an interval (see CsfGeneCoordsReader).
    """

    proteins: tuple[DecodeSomascanProtein, ...]
    conflict_distance_bp: int = SelectionConfig().cis_window_bp

    def __call__(self, table: pl.DataFrame) -> dict[ProteinId, GenomicInterval]:
        genes = magma_genes(table)
        intervals = {
            ProteinId(protein.seq_id): gene_interval(
                _gene_keys(protein), genes, self.conflict_distance_bp
            )
            for protein in self.proteins
        }
        return {pid: iv for pid, iv in intervals.items() if iv is not None}
