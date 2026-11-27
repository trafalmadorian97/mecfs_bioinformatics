"""
Source-specific readers turning a proteomics source's protein table into the UniProt accessions
of each measured protein: one row per distinct (protein id, UniProt accession).

A protein complex (for example EBI3_IL27) gives one row per subunit accession.
"""

from typing import Protocol

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.decode_database.decode_somascan_protein import (
    DecodeSomascanProtein,
)
from mecfs_bio.build_system.task.ppp_ldsc.gene_coords import (
    ST3_OLINK_ID_COL,
    ST3_UNIPROT_COL,
)

PROTEIN_ACCESSIONS_PROTEIN_ID_COL = "protein_id"
PROTEIN_ACCESSIONS_UNIPROT_COL = "uniprot_id"
_ST3_SEPARATOR = ";"
# ST29 lists a multi-protein target's accessions comma-separated.
_ST29_SEPARATOR = ","

PROTEIN_ACCESSIONS_SCHEMA = pl.Schema(
    {
        PROTEIN_ACCESSIONS_PROTEIN_ID_COL: pl.String,
        PROTEIN_ACCESSIONS_UNIPROT_COL: pl.String,
    }
)


class ProteinAccessionsReader(Protocol):
    def __call__(self, table: pl.DataFrame) -> pl.DataFrame:
        """One row per distinct (protein id, UniProt accession), in PROTEIN_ACCESSIONS_SCHEMA."""
        ...


@frozen(slots=True)
class St3ProteinAccessionsReader:
    """UKB-PPP: the extracted Sun et al. 2023 ST3 sheet, keyed by Olink id."""

    def __call__(self, table: pl.DataFrame) -> pl.DataFrame:
        return (
            table.select(
                pl.col(ST3_OLINK_ID_COL)
                .cast(pl.String)
                .alias(PROTEIN_ACCESSIONS_PROTEIN_ID_COL),
                pl.col(ST3_UNIPROT_COL)
                .cast(pl.String)
                .str.split(_ST3_SEPARATOR)
                .alias(PROTEIN_ACCESSIONS_UNIPROT_COL),
            )
            .explode(PROTEIN_ACCESSIONS_UNIPROT_COL)
            .with_columns(pl.col(PROTEIN_ACCESSIONS_UNIPROT_COL).str.strip_chars())
            .unique(maintain_order=True)
            .cast(PROTEIN_ACCESSIONS_SCHEMA)
        )


@frozen(slots=True)
class DecodeProteinAccessionsReader:
    """
    deCODE SomaScan: the UniProt accessions of the manifest rows it holds, keyed by seq id.
    The table it is called with is ignored.
    """

    proteins: tuple[DecodeSomascanProtein, ...]

    def __call__(self, table: pl.DataFrame) -> pl.DataFrame:
        return (
            pl.DataFrame(
                {
                    PROTEIN_ACCESSIONS_PROTEIN_ID_COL: [
                        p.seq_id for p in self.proteins
                    ],
                    PROTEIN_ACCESSIONS_UNIPROT_COL: [p.uniprot for p in self.proteins],
                },
                schema=PROTEIN_ACCESSIONS_SCHEMA,
            )
            .with_columns(
                pl.col(PROTEIN_ACCESSIONS_UNIPROT_COL).str.split(_ST29_SEPARATOR)
            )
            .explode(PROTEIN_ACCESSIONS_UNIPROT_COL)
            .with_columns(pl.col(PROTEIN_ACCESSIONS_UNIPROT_COL).str.strip_chars())
            .unique(maintain_order=True)
        )
