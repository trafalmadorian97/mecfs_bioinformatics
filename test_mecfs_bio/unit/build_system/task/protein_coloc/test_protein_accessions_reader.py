import polars as pl

from mecfs_bio.build_system.task.ppp_ldsc.gene_coords import (
    ST3_OLINK_ID_COL,
    ST3_UNIPROT_COL,
)
from mecfs_bio.build_system.task.protein_coloc.protein_accessions_reader import (
    PROTEIN_ACCESSIONS_SCHEMA,
    St3ProteinAccessionsReader,
)


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
