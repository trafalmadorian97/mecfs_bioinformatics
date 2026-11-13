"""
Source-specific readers turning a proteomics source's protein table into the genes behind each
measured protein: one row per (protein id, gene symbol, UniProt accession).

A protein complex (for example EBI3_IL27) gives one row per gene. Genes duplicated in the genome
that encode one protein (for example DEFA1 and DEFA1B, both P59665) keep a row per symbol.
"""

from typing import Protocol

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.ppp_ldsc.gene_coords import (
    ST3_GENE_SYMBOL_COL,
    ST3_OLINK_ID_COL,
    ST3_UNIPROT_COL,
)

PROTEIN_GENES_PROTEIN_ID_COL = "protein_id"
PROTEIN_GENES_SYMBOL_COL = "gene_symbol"
PROTEIN_GENES_UNIPROT_COL = "uniprot_id"
_ST3_SEPARATOR = ";"

PROTEIN_GENES_SCHEMA = pl.Schema(
    {
        PROTEIN_GENES_PROTEIN_ID_COL: pl.String,
        PROTEIN_GENES_SYMBOL_COL: pl.String,
        PROTEIN_GENES_UNIPROT_COL: pl.String,
    }
)


class ProteinGenesReader(Protocol):
    def __call__(self, table: pl.DataFrame) -> pl.DataFrame:
        """One row per (protein id, gene symbol, UniProt accession), in PROTEIN_GENES_SCHEMA."""
        ...


@frozen(slots=True)
class St3ProteinGenesReader:
    """UKB-PPP: the extracted Sun et al. 2023 ST3 sheet, keyed by Olink id."""

    def __call__(self, table: pl.DataFrame) -> pl.DataFrame:
        split = table.select(
            pl.col(ST3_OLINK_ID_COL)
            .cast(pl.String)
            .alias(PROTEIN_GENES_PROTEIN_ID_COL),
            pl.col(ST3_GENE_SYMBOL_COL)
            .str.split(_ST3_SEPARATOR)
            .alias(PROTEIN_GENES_SYMBOL_COL),
            pl.col(ST3_UNIPROT_COL)
            .str.split(_ST3_SEPARATOR)
            .alias(PROTEIN_GENES_UNIPROT_COL),
        )
        mismatched = split.filter(
            pl.col(PROTEIN_GENES_SYMBOL_COL).list.len()
            != pl.col(PROTEIN_GENES_UNIPROT_COL).list.len()
        )
        assert mismatched.height == 0, (
            f"gene symbol and UniProt lists differ in length for "
            f"{mismatched[PROTEIN_GENES_PROTEIN_ID_COL].to_list()}"
        )
        return (
            split.explode(PROTEIN_GENES_SYMBOL_COL, PROTEIN_GENES_UNIPROT_COL)
            .with_columns(
                pl.col(PROTEIN_GENES_SYMBOL_COL).str.strip_chars(),
                pl.col(PROTEIN_GENES_UNIPROT_COL).str.strip_chars(),
            )
            .unique(maintain_order=True)
            .cast(PROTEIN_GENES_SCHEMA)
        )
