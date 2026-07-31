"""
Check that reading ST3 gene coordinates from an extracted dataframe reproduces exactly what the
previous openpyxl reader produced from the spreadsheet.

The LDSC tasks moved from opening the Sun et al. supplementary xlsx directly to consuming the
sheet already extracted to parquet, so that they can be unit tested without a spreadsheet fixture.
That refactor is only safe if the two readers agree on every protein: a silently different set of
cis windows would change every cis-excluded heritability and genetic correlation without failing
anything. This script reimplements the old openpyxl reader verbatim and diffs it against the new
dataframe path.

Run: pixi r python experiments/claude/ppp_ldsc/verify_st3_gene_coords_refactor.py
"""

from pathlib import Path

import openpyxl
import pandas as pd
import polars as pl

from mecfs_bio.build_system.task.ppp_ldsc.gene_coords import (
    ST3_GENE_CHROM_COL,
    ST3_GENE_END_COL,
    ST3_GENE_START_COL,
    ST3_OLINK_ID_COL,
    parse_chrom,
    read_gene_coords,
)
from mecfs_bio.constants.ppp_database_constants import Oid
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

XLSX_PATH = Path(
    "assets/base_asset_store/reference_data/pqtl_data/sun_et_al_2023/raw/"
    "sun_et_al_2023_supplementary.xlsx"
)
_ST3_SHEET = "ST3"


def read_gene_coords_via_openpyxl(xlsx_path: Path) -> dict[Oid, GenomicInterval]:
    """The reader as it stood before the refactor, kept here as the reference."""
    workbook = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    sheet = workbook[_ST3_SHEET]
    header_index: dict[str, int] | None = None
    coords: dict[Oid, GenomicInterval] = {}
    for row in sheet.iter_rows(values_only=True):
        if header_index is None:
            if row and ST3_OLINK_ID_COL in row:
                header_index = {
                    str(name): i for i, name in enumerate(row) if name is not None
                }
            continue
        oid = row[header_index[ST3_OLINK_ID_COL]]
        chrom = parse_chrom(row[header_index[ST3_GENE_CHROM_COL]])
        start = row[header_index[ST3_GENE_START_COL]]
        end = row[header_index[ST3_GENE_END_COL]]
        if oid is None or chrom is None or start is None or end is None:
            continue
        coords[Oid(oid)] = GenomicInterval(chrom=chrom, start=int(start), end=int(end))
    assert header_index is not None
    return coords


def extract_st3_as_dataframe(xlsx_path: Path) -> pl.DataFrame:
    """The same extraction ExtractSheetFromExelFileTask performs when wired for ST3."""
    frame = pd.read_excel(xlsx_path, sheet_name=_ST3_SHEET, skiprows=2)
    for column in (ST3_GENE_CHROM_COL, ST3_GENE_START_COL, ST3_GENE_END_COL):
        frame[column] = frame[column].astype(str)
    return pl.from_pandas(frame)


def main() -> None:
    reference = read_gene_coords_via_openpyxl(XLSX_PATH)
    refactored = read_gene_coords(extract_st3_as_dataframe(XLSX_PATH))

    print(f"openpyxl reader:  {len(reference)} proteins with coordinates")
    print(f"dataframe reader: {len(refactored)} proteins with coordinates")

    only_reference = sorted(set(reference) - set(refactored))
    only_refactored = sorted(set(refactored) - set(reference))
    differing = sorted(
        oid
        for oid in set(reference) & set(refactored)
        if reference[oid] != refactored[oid]
    )
    print(f"only in openpyxl reader:  {only_reference[:10]} (n={len(only_reference)})")
    print(f"only in dataframe reader: {only_refactored[:10]} (n={len(only_refactored)})")
    print(f"differing coordinates:    {differing[:10]} (n={len(differing)})")

    assert reference == refactored, "the two readers disagree; see the diffs printed above"
    print("MATCH: both readers produce identical gene coordinates")


if __name__ == "__main__":
    main()
