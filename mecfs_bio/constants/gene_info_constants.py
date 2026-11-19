"""
Columns that AnnotateGenesFromReferenceTask adds to a table of Ensembl gene ids.

A gene can have several values of a column (an Ensembl id that HGNC lists under two symbols, or a
gene encoding several UniProt proteins); they are joined by GENE_INFO_SEPARATOR.
"""

from typing import Final

GENE_INFO_SYMBOL_COL: Final = "gene_symbol"
GENE_INFO_NAME_COL: Final = "gene_name"
GENE_INFO_UNIPROT_ID_COL: Final = "uniprot_id"
GENE_INFO_PROTEIN_NAMES_COL: Final = "protein_names"
GENE_INFO_UNIPROT_DESCRIPTION_COL: Final = "uniprot_description"
GENE_INFO_SUBCELLULAR_LOCALISATION_COL: Final = "subcellular_localisation"

GENE_INFO_COLS: Final = (
    GENE_INFO_SYMBOL_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_UNIPROT_ID_COL,
    GENE_INFO_PROTEIN_NAMES_COL,
    GENE_INFO_UNIPROT_DESCRIPTION_COL,
    GENE_INFO_SUBCELLULAR_LOCALISATION_COL,
)

GENE_INFO_SEPARATOR: Final = " | "
