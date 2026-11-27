"""
One deCODE SomaScan assay's summary statistics file (a manifest row).
"""

from attrs import frozen

from mecfs_bio.build_system.wf.s3_access import S3ObjectRef
from mecfs_bio.constants.csf_database_constants import EntrezGeneId, SeqId, UniProtId


@frozen(slots=True)
class DecodeSomascanProtein:
    """
    seq_id: the SomaScan sequence id in dash form (e.g. 10076-1), the assay's primary key.
    gene_symbol, target_name: from the file name; the symbol labels assets but is not unique
        across assays and may be an outdated HGNC symbol.
    uniprot: the target's UniProt accession, from ST29.
    entrez_gene_id: NCBI Entrez gene id(s) of the accession via HGNC, "|"-joined when several;
        None when HGNC gives none.
    ref: the file's S3 object; size and md5 are its expected size in bytes and md5.
    """

    seq_id: SeqId
    gene_symbol: str
    target_name: str
    uniprot: UniProtId
    entrez_gene_id: EntrezGeneId | None
    ref: S3ObjectRef
    size: int
    md5: str

    @property
    def key(self) -> str:
        return self.ref.key
