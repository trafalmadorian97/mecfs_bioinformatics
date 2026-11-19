"""
The reference tables used to annotate gene lists with gene symbols, names and UniProt protein
information.
"""

from mecfs_bio.assets.reference_data.hgnc.hgnc_complete_set import (
    HGNC_ENSEMBL_GENES,
    HGNC_UNIPROT_GENES,
)
from mecfs_bio.assets.reference_data.uniprot.uniprot_lookup_table import UNIPROT_LOOKUP
from mecfs_bio.build_system.task.annotate_genes_from_reference_task import (
    GeneInfoReferenceTasks,
)

GENE_INFO_REFERENCE = GeneInfoReferenceTasks(
    hgnc_genes_task=HGNC_ENSEMBL_GENES,
    uniprot_genes_task=HGNC_UNIPROT_GENES,
    uniprot_task=UNIPROT_LOOKUP,
)
