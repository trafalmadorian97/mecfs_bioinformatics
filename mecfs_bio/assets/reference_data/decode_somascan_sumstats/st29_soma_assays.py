"""
Read the SomaScan assays of ST29 from the supplementary tables of Eldjarn et al. 2023,
"Large-scale plasma proteomics comparisons through genetics and disease associations",
Nature 622, 348-358, doi 10.1038/s41586-023-06563-x, for regenerate_decode_manifest.py.

Source: sheet ST29_protein_classification of the workbook (title on row 1, header on row 3).
ST29 has one row per protein targeted by Olink or SomaScan; the rows with a SeqId and
"Targeted Soma" == "Y" are the 4,907 SomaScan assays with deCODE summary statistics.

ST29 repeats a SomaScan assay once per Olink assay of the same protein (4,930 rows for 4,907
seq ids), and its confidence tier and gene label belong to the pairing. st29_soma_assays keeps
one row per seq id: the best (lowest) confidence tier, and the gene of the best-tier row (ties
broken alphabetically, which picks NPPB over NTproBNP). UniProt and cis pQTL soma are
per-assay and are asserted to agree across the repeats. Its columns are seq_id (dash form, e.g.
10001-7; ST29 writes 10001_7), uniprot (as ST29 writes it, several accessions comma
separated), gene, cis_pqtl_soma and confidence_tier.
"""

from pathlib import Path

import openpyxl
import polars as pl

_ST29_SHEET = "ST29_protein_classification"
_HEADER_ROW = 3
EXPECTED_SOMA_ASSAYS = 4_907

ST29_SEQ_ID_COL = "seq_id"
ST29_UNIPROT_COL = "uniprot"
ST29_GENE_COL = "gene"
ST29_CIS_PQTL_SOMA_COL = "cis_pqtl_soma"
ST29_CONFIDENCE_TIER_COL = "confidence_tier"

_SOURCE_COLUMNS = {
    "SeqId": ST29_SEQ_ID_COL,
    "UniProt": ST29_UNIPROT_COL,
    "Gene": ST29_GENE_COL,
    "cis pQTL soma": ST29_CIS_PQTL_SOMA_COL,
    "confidence tier": ST29_CONFIDENCE_TIER_COL,
}


def read_sheet(workbook: Path, sheet: str, header_row: int) -> pl.DataFrame:
    """All rows below header_row of one sheet, every cell as a string (None when empty)."""
    wb = openpyxl.load_workbook(workbook, read_only=True)
    rows = list(wb[sheet].iter_rows(min_row=header_row, values_only=True))
    header = [str(h) for h in rows[0]]
    body = [[None if v is None else str(v) for v in row] for row in rows[1:]]
    return pl.DataFrame(body, schema=dict.fromkeys(header, pl.String), orient="row")


def soma_assays(st29: pl.DataFrame) -> pl.DataFrame:
    return (
        st29.filter(pl.col("SeqId").is_not_null() & (pl.col("Targeted Soma") == "Y"))
        .select(list(_SOURCE_COLUMNS))
        .rename(_SOURCE_COLUMNS)
        .with_columns(
            pl.col(ST29_SEQ_ID_COL).str.replace("_", "-", literal=True),
            pl.col(ST29_CONFIDENCE_TIER_COL).cast(pl.Int32),
        )
    )


def one_row_per_assay(assays: pl.DataFrame) -> pl.DataFrame:
    per_assay = [ST29_UNIPROT_COL, ST29_CIS_PQTL_SOMA_COL]
    conflicts = (
        assays.group_by(ST29_SEQ_ID_COL)
        .agg(pl.col(c).n_unique().alias(c) for c in per_assay)
        .filter(pl.any_horizontal(pl.col(c) > 1 for c in per_assay))
    )
    assert conflicts.is_empty(), f"per-assay ST29 columns disagree:\n{conflicts}"
    return (
        assays.sort(ST29_SEQ_ID_COL, ST29_CONFIDENCE_TIER_COL, ST29_GENE_COL)
        .group_by(ST29_SEQ_ID_COL, maintain_order=True)
        .first()
        .select(list(_SOURCE_COLUMNS.values()))
    )


def st29_soma_assays(workbook: Path) -> pl.DataFrame:
    """One row per SomaScan assay of ST29, read from the supplementary workbook."""
    assays = one_row_per_assay(
        soma_assays(read_sheet(workbook, _ST29_SHEET, _HEADER_ROW))
    )
    assert assays.height == EXPECTED_SOMA_ASSAYS, assays.height
    assert assays[ST29_SEQ_ID_COL].is_unique().all()
    return assays
