import polars as pl
import pytest

from mecfs_bio.build_system.task.ppp_ldsc.gene_coords import (
    ST3_GENE_SYMBOL_COL,
    ST3_OLINK_ID_COL,
    ST3_UNIPROT_COL,
)
from mecfs_bio.build_system.task.protein_coloc.protein_genes_reader import (
    PROTEIN_GENES_SCHEMA,
    St3ProteinGenesReader,
)


def _st3(rows: list[tuple[str, str, str]]) -> pl.DataFrame:
    return pl.DataFrame(
        {
            ST3_OLINK_ID_COL: [r[0] for r in rows],
            ST3_GENE_SYMBOL_COL: [r[1] for r in rows],
            ST3_UNIPROT_COL: [r[2] for r in rows],
        }
    )


def test_complex_gives_one_row_per_gene() -> None:
    genes = St3ProteinGenesReader()(
        _st3(
            [
                ("OID1", "AAA", "P00001"),
                ("OID2", "EBI3;IL27", "Q14213;Q8NEV9"),
                ("OID3", "DEFA1;DEFA1B", "P59665;P59665"),
            ]
        )
    )
    assert genes.schema == PROTEIN_GENES_SCHEMA
    assert genes.rows() == [
        ("OID1", "AAA", "P00001"),
        ("OID2", "EBI3", "Q14213"),
        ("OID2", "IL27", "Q8NEV9"),
        ("OID3", "DEFA1", "P59665"),
        ("OID3", "DEFA1B", "P59665"),
    ]


def test_mismatched_list_lengths_fail() -> None:
    with pytest.raises(AssertionError):
        St3ProteinGenesReader()(_st3([("OID1", "AAA;BBB", "P00001")]))
