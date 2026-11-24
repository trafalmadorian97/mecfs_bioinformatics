"""
Extract the SomaScan rows of ST29 from the supplementary tables of Eldjarn et al. 2023,
"Large-scale plasma proteomics comparisons through genetics and disease associations",
Nature 622, 348-358, doi 10.1038/s41586-023-06563-x.

Source: Supplementary Tables, file 41586_2023_6563_MOESM3_ESM.xlsx, sheet
ST29_protein_classification (title on row 1, header on row 3). ST29 has one row per protein
targeted by Olink or SomaScan; the rows with a SeqId and "Targeted Soma" == "Y" are the 4,907
SomaScan assays with deCODE summary statistics.

ST29 repeats a SomaScan assay once per Olink assay of the same protein (4,930 rows for 4,907
seq ids), and its confidence tier and gene label belong to the pairing. The extract keeps one
row per seq id: the best (lowest) confidence tier, and the gene of the best-tier row (ties
broken alphabetically, which picks NPPB over NTproBNP). UniProt and cis pQTL soma are
per-assay and are asserted to agree across the repeats.

Writes st29_soma_assays.csv next to this script, with columns seq_id (dash form, e.g.
10001-7; ST29 writes 10001_7), uniprot (as ST29 writes it, several accessions comma
separated), gene, cis_pqtl_soma and confidence_tier.

With --st36, extracts sheet ST36_pleiotropic_pQTLs instead: one row per highly pleiotropic
pQTL (Olink, SomaScan SMP or SomaScan PC0), with the chromosome as a code (chrX = 23), the
Olink, SMP and PC0 lead positions (hg38), the three pleiotropy flags as booleans, and the SMP
lead's name, rsName and number of proteins. Writes st36_pleiotropic_pqtls.csv.

Usage (from the repo root), with the workbook downloaded from the paper:
  pixi r python mecfs_bio/assets/reference_data/decode_somascan_sumstats/extract_st29.py \
      [--st36] <41586_2023_6563_MOESM3_ESM.xlsx>
"""

import sys
from pathlib import Path

import openpyxl
import polars as pl

from mecfs_bio.constants.decode_somascan_constants import (
    DECODE_CHR_X_CODE,
    ST36_CHROM_COL,
    ST36_OLINK_POS_COL,
    ST36_PC0_POS_COL,
    ST36_PLEIOTROPIC_OLINK_COL,
    ST36_PLEIOTROPIC_PC0_COL,
    ST36_PLEIOTROPIC_SMP_COL,
    ST36_SMP_NAME_COL,
    ST36_SMP_NPROT_COL,
    ST36_SMP_POS_COL,
    ST36_SMP_RSNAME_COL,
)

_DIR = Path(__file__).parent
ST29_CSV_PATH = _DIR / "st29_soma_assays.csv"
_ST29_SHEET = "ST29_protein_classification"
_HEADER_ROW = 3
ST36_CSV_PATH = _DIR / "st36_pleiotropic_pqtls.csv"
_ST36_SHEET = "ST36_pleiotropic_pQTLs"
# ST36 headers wrap onto several lines; they are matched with newlines replaced by spaces.
_ST36_COLUMNS = {
    "Chrom": ST36_CHROM_COL,
    "Olink Pos (hg38)": ST36_OLINK_POS_COL,
    "Highly Pleiotropic Olink": ST36_PLEIOTROPIC_OLINK_COL,
    "Highly Pleiotropic Soma SMP": ST36_PLEIOTROPIC_SMP_COL,
    "Highly Pleiotropic Soma PC0": ST36_PLEIOTROPIC_PC0_COL,
    "Soma SMP Pos": ST36_SMP_POS_COL,
    "Soma SMP Name": ST36_SMP_NAME_COL,
    "Soma SMP rsName": ST36_SMP_RSNAME_COL,
    "Soma SMP nprot": ST36_SMP_NPROT_COL,
    "Soma PC0 Pos": ST36_PC0_POS_COL,
}
_ST36_FLAGS = (
    ST36_PLEIOTROPIC_OLINK_COL,
    ST36_PLEIOTROPIC_SMP_COL,
    ST36_PLEIOTROPIC_PC0_COL,
)
_ST36_SOURCE_FLAGS = [k for k, v in _ST36_COLUMNS.items() if v in _ST36_FLAGS]
_ST36_POSITIONS = (ST36_OLINK_POS_COL, ST36_SMP_POS_COL, ST36_PC0_POS_COL)
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


def st36_pleiotropic_pqtls(st36: pl.DataFrame) -> pl.DataFrame:
    renamed = st36.rename(lambda c: c.replace("\n", " ")).select(list(_ST36_COLUMNS))
    flag_values = set(
        renamed.select(pl.concat_list(_ST36_SOURCE_FLAGS).explode())
        .to_series()
        .drop_nulls()
        .to_list()
    )
    assert flag_values <= {"Y", "N"}, f"unexpected ST36 flag values {flag_values}"
    table = renamed.rename(_ST36_COLUMNS).with_columns(
        pl.col(ST36_CHROM_COL)
        .str.strip_prefix("chr")
        .replace("X", str(DECODE_CHR_X_CODE))
        .cast(pl.Int32),
        *(pl.col(c).cast(pl.Float64).cast(pl.Int64) for c in _ST36_POSITIONS),
        pl.col(ST36_SMP_NPROT_COL).cast(pl.Float64).cast(pl.Int32),
        *((pl.col(c) == "Y").fill_null(False) for c in _ST36_FLAGS),
    )
    smp = table.filter(pl.col(ST36_PLEIOTROPIC_SMP_COL))
    assert smp[ST36_SMP_POS_COL].null_count() == 0, (
        "an SMP-pleiotropic row has no SMP Pos"
    )
    return table


def main() -> None:
    if sys.argv[1] == "--st36":
        st36 = st36_pleiotropic_pqtls(
            read_sheet(Path(sys.argv[2]), _ST36_SHEET, _HEADER_ROW)
        )
        st36.write_csv(ST36_CSV_PATH)
        print(f"wrote {st36.height} rows to {ST36_CSV_PATH}")
        return
    workbook = Path(sys.argv[1])
    assays = one_row_per_assay(
        soma_assays(read_sheet(workbook, _ST29_SHEET, _HEADER_ROW))
    )
    assert assays.height == EXPECTED_SOMA_ASSAYS, assays.height
    assert assays[ST29_SEQ_ID_COL].is_unique().all()
    assays.write_csv(ST29_CSV_PATH)
    print(f"wrote {assays.height} rows to {ST29_CSV_PATH}")


if __name__ == "__main__":
    main()
