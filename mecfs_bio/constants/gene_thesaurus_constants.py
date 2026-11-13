"""
Column names of the gene thesaurus (an Ensembl BioMart export of all human genes: Ensembl gene
ids, gene symbols, descriptions and coordinates).
"""

from typing import Final

THESAURUS_ENSEMBL_ID_COL: Final = "Gene stable ID"
THESAURUS_GENE_NAME_COL: Final = "Gene name"
THESAURUS_DESCRIPTION_COL: Final = "Gene description"
THESAURUS_CHROM_COL: Final = "Chromosome/scaffold name"

# Chromosomes of the primary assembly. Genes on alternate-haplotype or patch scaffolds (for
# example HSCHR2_1_CTG5) duplicate primary-assembly genes under the same symbol.
THESAURUS_PRIMARY_CHROMOSOMES: Final = frozenset(
    [str(c) for c in range(1, 23)] + ["X", "Y", "MT"]
)
