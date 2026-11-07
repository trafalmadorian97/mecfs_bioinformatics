"""Per-trait summary of CSF coloc raw results: loci, candidates, status counts, failure reasons,
and the pairs with PP.H4 >= 0.8.
Usage: pixi r python experiments/claude/csf_coloc/coloc_run_summary.py <name>..."""

import sys

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.analysis import csf_coloc_analysis as a
from mecfs_bio.asset_generator.protein_coloc_asset_generator import ColocTaskGroup
from mecfs_bio.constants.protein_coloc_constants import (
    CIS_STATUS_COL,
    FAILURE_REASON_COL,
    GENE_COL,
    LOCUS_ID_COL,
    PP_H3_COL,
    PP_H4_COL,
    PROTEIN_ID_COL,
    RAW_RESULTS_FILENAME,
    STATUS_COL,
)

H4 = 0.8


GROUPS: dict[str, ColocTaskGroup] = {
    "ad": a.BELLENGUEZ_ALZHEIMERS_CSF_COLOC,
    "decode_me": a.DECODE_ME_CSF_COLOC,
    "ra": a.RA_SEROPOSITIVE_CSF_COLOC,
    "ldl": a.WILLER_LDL_CSF_COLOC,
    "ibd": a.LIU_ET_AL_2023_IBD_CSF_COLOC,
    "scz": a.PGC_2022_SCH_CSF_COLOC,
    "pd": a.NALLS_2019_PARKINSONS_CSF_COLOC,
    "als": a.VAN_RHEENEN_2021_ALS_CSF_COLOC,
    "decode_me_1e6": a.DECODE_ME_CSF_COLOC_LEAD_P_1E_6,
    "decode_me_1e5": a.DECODE_ME_CSF_COLOC_LEAD_P_1E_5,
}


def summarize(name: str, group: ColocTaskGroup) -> None:
    out_dir = DEFAULT_RUNNER.meta_to_path(group.coloc_task.meta)
    raw = pl.read_parquet(out_dir / RAW_RESULTS_FILENAME)
    print(f"===== {name} ({group.coloc_task.asset_id})")
    print(
        f"loci with candidates: {raw[LOCUS_ID_COL].n_unique()}; result rows: {raw.height}"
    )
    print(raw.group_by(STATUS_COL).len().sort(STATUS_COL))
    failed = raw.filter(pl.col(FAILURE_REASON_COL).is_not_null())
    if failed.height:
        print(failed.group_by(FAILURE_REASON_COL).len())
    hits = raw.filter(pl.col(PP_H4_COL) >= H4).sort(PP_H4_COL, descending=True)
    print(f"PP.H4 >= {H4}: {hits.height} pairs")
    print(
        hits.select(
            LOCUS_ID_COL, PROTEIN_ID_COL, GENE_COL, CIS_STATUS_COL, PP_H3_COL, PP_H4_COL
        )
    )


def main() -> None:
    pl.Config.set_tbl_rows(60)
    pl.Config.set_tbl_width_chars(200)
    for name in [n for n in sys.argv[1:] if n != "--top"]:
        summarize(name, GROUPS[name])
        if "--top" in sys.argv:
            top_pairs(name)


def top_pairs(name: str, n: int = 10) -> None:
    """The n tested pairs with the highest PP.H4, and every pair in the extended MHC."""
    out_dir = DEFAULT_RUNNER.meta_to_path(GROUPS[name].coloc_task.meta)
    raw = pl.read_parquet(out_dir / RAW_RESULTS_FILENAME)
    cols = [
        LOCUS_ID_COL,
        PROTEIN_ID_COL,
        GENE_COL,
        CIS_STATUS_COL,
        STATUS_COL,
        PP_H3_COL,
        PP_H4_COL,
        "LEAD_MLOG10P",
    ]
    print(f"===== {name}: top {n} by PP.H4")
    print(raw.sort(PP_H4_COL, descending=True, nulls_last=True).head(n).select(cols))
    print(f"===== {name}: chr6 25-34 Mb")
    print(
        raw.filter(
            pl.col("CHR") == 6,
            pl.col("END") >= 25_000_000,
            pl.col("START") <= 34_000_000,
        ).select(cols)
    )


if __name__ == "__main__":
    main()
