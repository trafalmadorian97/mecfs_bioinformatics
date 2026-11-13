"""
Column names and values of the HGNC complete set (https://www.genenames.org/download/): one row
per approved human gene, with its current symbol, its Ensembl gene id and the UniProt accessions
of its protein products.
"""

from typing import Final

HGNC_SYMBOL_COL: Final = "symbol"
HGNC_NAME_COL: Final = "name"
HGNC_LOCUS_TYPE_COL: Final = "locus_type"
HGNC_ENSEMBL_ID_COL: Final = "ensembl_gene_id"
# "|"-separated UniProt accessions.
HGNC_UNIPROT_IDS_COL: Final = "uniprot_ids"
HGNC_LIST_SEPARATOR: Final = "|"

# A transcript spanning two adjacent genes (for example PTGES3L-AARSD1). HGNC lists the
# readthrough under the UniProt accession of one of the genes it spans.
HGNC_READTHROUGH_LOCUS_TYPE: Final = "readthrough"
