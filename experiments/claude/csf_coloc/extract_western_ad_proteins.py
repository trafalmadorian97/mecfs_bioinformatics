"""
Extract Western et al. 2024's AD colocalization reference list from their Supplementary Tables
(Nature Genetics 56:2672, PMID 39528825; file 41588_2024_1972_MOESM4_ESM.xlsx from
static-content.springer.com) into western_ad_proteins.csv:

- ST43 (single-variant coloc.abf of each index pQTL with the AD GWAS, coloc v5.2.2): every
  aptamer x index pQTL pair with PP.H4 >= 0.8 (source "ST43").
- ST45 (aptamers associated with AD by at least two of PWAS, COLOC, MR): every row (source
  "ST45"), with their PP.H4 and which methods agree.

Usage: pixi r python experiments/claude/csf_coloc/extract_western_ad_proteins.py <MOESM4.xlsx>
"""

import sys
from pathlib import Path

import openpyxl
import polars as pl

OUT = Path(__file__).parent / "western_ad_proteins.csv"


def _table(path: str, sheet: str, ncol: int) -> pl.DataFrame:
    wb = openpyxl.load_workbook(path, read_only=True)
    rows = list(wb[sheet].iter_rows(values_only=True))
    header = [str(h) for h in rows[1][:ncol]]
    records = [dict(zip(header, r[:ncol])) for r in rows[2:] if r[0]]
    return pl.DataFrame(records, infer_schema_length=None).cast(pl.String)


def _seq_id(analyte: pl.Expr) -> pl.Expr:
    return analyte.str.strip_prefix("X").str.replace(".", "-", literal=True)


def main() -> None:
    path = sys.argv[1]
    st43 = (
        _table(path, "Supplementary Table 43", 11)
        .filter(pl.col("PP.H4.abf").cast(pl.Float64) >= 0.8)
        .select(
            pl.lit("ST43").alias("source"),
            pl.col("SomaID").alias("analyte"),
            pl.col("EntrezGeneSymbol").alias("gene"),
            pl.col("SNPID").alias("pqtl_id"),
            pl.col("PP.H4.abf").cast(pl.Float64).alias("their_pp_h4"),
            pl.col("pQTL Type").alias("pqtl_type"),
            pl.lit(None, dtype=pl.String).alias("methods"),
        )
    )
    st45 = _table(path, "Supplementary Table 45", 12).select(
        pl.lit("ST45").alias("source"),
        pl.col("SomaID").alias("analyte"),
        pl.col("EntrezGeneSymbol").alias("gene"),
        pl.col("PQTL.ID").alias("pqtl_id"),
        pl.col("PP.H4").cast(pl.Float64).alias("their_pp_h4"),
        pl.col("pQTL Type").alias("pqtl_type"),
        pl.col("Overlap").alias("methods"),
    )
    out = pl.concat([st43, st45]).with_columns(
        _seq_id(pl.col("analyte")).alias("seq_id")
    )
    out.write_csv(OUT)
    print(f"wrote {out.height} rows ({st43.height} ST43, {st45.height} ST45) to {OUT}")


if __name__ == "__main__":
    main()
