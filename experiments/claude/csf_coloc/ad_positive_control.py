"""
Positive-control check of the Alzheimer's x CSF coloc run against Western et al. 2024's own AD
colocalization (western_ad_proteins.csv, from their Supplementary Tables 43 and 45).

For every aptamer in their list, report from our raw results: whether any of its pairs was a
candidate (its hit segment overlapped an AD locus window), whether it was selected within the
pair budget, whether it was tested, and our best PP.H4; plus whether our locus window is in the
excluded MHC. Then list our PP.H4 >= 0.8 pairs that are not in their list.

Usage: pixi r python experiments/claude/csf_coloc/ad_positive_control.py
"""

from pathlib import Path

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.analysis.bellenguez_alzheimers_csf_coloc import (
    BELLENGUEZ_ALZHEIMERS_CSF_COLOC,
)
from mecfs_bio.constants.protein_coloc_constants import (
    CIS_STATUS_COL,
    GENE_COL,
    LOCUS_ID_COL,
    PP_H4_COL,
    PROTEIN_ID_COL,
    RAW_RESULTS_FILENAME,
    STATUS_COL,
    STATUS_FAILED,
    STATUS_SELECTED,
    STATUS_TESTED,
)

REFERENCE = Path(__file__).parent / "western_ad_proteins.csv"
H4 = 0.8


def _raw_results() -> pl.DataFrame:
    out_dir = DEFAULT_RUNNER.meta_to_path(
        BELLENGUEZ_ALZHEIMERS_CSF_COLOC.coloc_task.meta
    )
    return pl.read_parquet(out_dir / RAW_RESULTS_FILENAME)


def _per_aptamer(results: pl.DataFrame) -> pl.DataFrame:
    pairs = results.filter(pl.col(PROTEIN_ID_COL).is_not_null())
    selected = [STATUS_SELECTED, STATUS_TESTED, STATUS_FAILED]
    return pairs.group_by(PROTEIN_ID_COL).agg(
        pl.len().alias("our_candidate_pairs"),
        pl.col(STATUS_COL).is_in(selected).sum().alias("our_selected_pairs"),
        (pl.col(STATUS_COL) == STATUS_TESTED).sum().alias("our_tested_pairs"),
        (pl.col(STATUS_COL) == STATUS_FAILED).sum().alias("our_failed_pairs"),
        pl.col(PP_H4_COL).max().alias("our_max_pp_h4"),
        pl.col(LOCUS_ID_COL)
        .sort_by(PP_H4_COL, descending=True, nulls_last=True)
        .first()
        .alias("our_best_locus"),
        pl.col(CIS_STATUS_COL)
        .sort_by(PP_H4_COL, descending=True, nulls_last=True)
        .first()
        .alias("our_cis_status"),
    )


def _verdict(row: dict) -> str:
    if row["our_candidate_pairs"] is None:
        return "not a candidate"
    if row["our_tested_pairs"] == 0:
        return (
            "candidate, not tested"
            if row["our_selected_pairs"] == 0
            else "selected, failed"
        )
    return "colocalized" if (row["our_max_pp_h4"] or 0) >= H4 else "tested, PP.H4 < 0.8"


def main() -> None:
    results = _raw_results()
    print("our raw results by status:", dict(results.group_by(STATUS_COL).len().rows()))
    print("our loci:", results[LOCUS_ID_COL].n_unique())
    reference = pl.read_csv(REFERENCE)
    per_aptamer = _per_aptamer(results)
    joined = (
        reference.group_by("seq_id", "gene")
        .agg(
            pl.col("their_pp_h4").max(),
            pl.col("source").unique().sort().str.join("+"),
            pl.col("pqtl_type").unique().sort().str.join("+"),
            pl.col("pqtl_id").unique().sort().str.join(";"),
        )
        .join(per_aptamer, left_on="seq_id", right_on=PROTEIN_ID_COL, how="left")
    )
    joined = joined.with_columns(
        pl.Series("verdict", [_verdict(r) for r in joined.iter_rows(named=True)])
    ).sort("verdict", "their_pp_h4", descending=[False, True])
    with pl.Config(tbl_rows=200, tbl_cols=20, fmt_str_lengths=40, tbl_width_chars=250):
        print(joined.drop("pqtl_id"))
        print(joined.group_by("verdict", "pqtl_type").len().sort("verdict"))
        st45 = joined.filter(pl.col("source").str.contains("ST45"))
        print(
            "ST45 (>=2 methods) by verdict:",
            dict(st45.group_by("verdict").len().rows()),
        )
        theirs = set(reference["seq_id"].to_list())
        extra = results.filter(
            pl.col(PP_H4_COL) >= H4, ~pl.col(PROTEIN_ID_COL).is_in(theirs)
        ).select(PROTEIN_ID_COL, GENE_COL, LOCUS_ID_COL, CIS_STATUS_COL, PP_H4_COL)
        print(f"our PP.H4 >= {H4} pairs not in their list: {extra.height}")
        print(extra.sort(PP_H4_COL, descending=True))


if __name__ == "__main__":
    main()
