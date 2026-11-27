import polars as pl

from mecfs_bio.build_system.task.decode_database.decode_somascan_protein import (
    DecodeSomascanProtein,
)
from mecfs_bio.build_system.task.ppp_ldsc.gene_coords import (
    ST3_OLINK_ID_COL,
    ST3_UNIPROT_COL,
)
from mecfs_bio.build_system.task.protein_coloc.protein_accessions_reader import (
    PROTEIN_ACCESSIONS_SCHEMA,
    DecodeProteinAccessionsReader,
    St3ProteinAccessionsReader,
)
from mecfs_bio.build_system.wf.s3_access import S3Endpoint, S3ObjectRef
from mecfs_bio.constants.csf_database_constants import SeqId, UniProtId


def test_one_row_per_distinct_accession_of_each_protein() -> None:
    accessions = St3ProteinAccessionsReader()(
        pl.DataFrame(
            {
                ST3_OLINK_ID_COL: ["OID1", "OID2", "OID3"],
                ST3_UNIPROT_COL: ["P00001", "Q14213;Q8NEV9", "P59665;P59665"],
            }
        )
    )
    assert accessions.schema == PROTEIN_ACCESSIONS_SCHEMA
    assert accessions.rows() == [
        ("OID1", "P00001"),
        ("OID2", "Q14213"),
        ("OID2", "Q8NEV9"),
        ("OID3", "P59665"),
    ]


def test_decode_accessions_split_multi_uniprot() -> None:
    endpoint = S3Endpoint(
        url="https://fake.example", region="r", credentials_profile="p"
    )
    proteins = tuple(
        DecodeSomascanProtein(
            seq_id=SeqId(seq_id),
            gene_symbol="G",
            target_name="G",
            uniprot=UniProtId(uniprot),
            entrez_gene_id=None,
            ref=S3ObjectRef(endpoint=endpoint, bucket="b", key=seq_id),
            size=1,
            md5="0" * 32,
        )
        for seq_id, uniprot in [("1-1", "P1, P2"), ("2-1", "P3")]
    )
    # The table passed in is ignored: the manifest rows carry the accessions.
    accessions = DecodeProteinAccessionsReader(proteins=proteins)(pl.DataFrame())
    assert accessions.schema == PROTEIN_ACCESSIONS_SCHEMA
    assert accessions.rows() == [("1-1", "P1"), ("1-1", "P2"), ("2-1", "P3")]
