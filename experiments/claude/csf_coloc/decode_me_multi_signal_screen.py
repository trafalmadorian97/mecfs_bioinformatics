"""Screen the tested DecodeME x CSF pairs for signs of more than one signal per window, i.e.
pairs where SuSiE-coloc might find a colocalization that coloc.abf (one causal variant per
trait) misses.

Uses marginal statistics only (no LD), from each pair's aligned window:
- TRAIT_AT_PROT_LEAD: trait -log10 p at the protein's lead variant, and the maximum within
  +/-25 kb of it (TRAIT_NEAR_PROT_LEAD). A real secondary trait signal under the pQTL is the
  case SuSiE-coloc rescues.
- PROT_AT_TRAIT_LEAD: protein -log10 p at the trait's lead variant (the mirror case).
- TRAIT_PEAKS / PROT_PEAKS: distance-clumped peaks (greedy, 250 kb, trait p < 1e-4, protein
  p < 1e-5). Without LD this over-counts in long LD blocks (e.g. the MHC); it is a hint only.

Usage: pixi r python experiments/claude/csf_coloc/decode_me_multi_signal_screen.py
"""

import numpy as np
import polars as pl

from mecfs_bio.analysis import csf_coloc_analysis as a
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    two_sided_mlog10p,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_POS_COL
from mecfs_bio.constants.protein_coloc_constants import (
    BETA_PROTEIN_COL,
    BETA_TRAIT_COL,
    CIS_STATUS_COL,
    GENE_COL,
    LOCUS_ID_COL,
    PP_H3_COL,
    PP_H4_COL,
    PROTEIN_ID_COL,
    RAW_RESULTS_FILENAME,
    SE_PROTEIN_COL,
    SE_TRAIT_COL,
    STATUS_COL,
    STATUS_TESTED,
    WINDOWS_DIRNAME,
    pair_file_stem,
)

RUNS = {
    "5e-8": a.DECODE_ME_CSF_COLOC,
    "1e-6": a.DECODE_ME_CSF_COLOC_LEAD_P_1E_6,
    "1e-5": a.DECODE_ME_CSF_COLOC_LEAD_P_1E_5,
}
CLUMP_KB = 250
NEAR_KB = 25


def _clump(pos: np.ndarray, mlog10p: np.ndarray, threshold: float) -> int:
    order = np.argsort(-np.nan_to_num(mlog10p, nan=-1.0))
    peaks: list[int] = []
    for i in order:
        if not mlog10p[i] >= threshold:
            break
        if all(abs(pos[i] - p) > CLUMP_KB * 1000 for p in peaks):
            peaks.append(int(pos[i]))
    return len(peaks)


def _pair_stats(window: pl.DataFrame) -> dict[str, float | int]:
    pos = window[GWASLAB_POS_COL].to_numpy()
    with np.errstate(invalid="ignore", divide="ignore"):
        t = two_sided_mlog10p((window[BETA_TRAIT_COL] / window[SE_TRAIT_COL]).to_numpy())
        p = two_sided_mlog10p(
            (window[BETA_PROTEIN_COL] / window[SE_PROTEIN_COL]).to_numpy()
        )
    t_lead, p_lead = int(np.nanargmax(t)), int(np.nanargmax(p))
    near = np.abs(pos - pos[p_lead]) <= NEAR_KB * 1000
    return {
        "TRAIT_MAX": float(t[t_lead]),
        "PROT_MAX": float(p[p_lead]),
        "LEAD_DIST_KB": abs(int(pos[t_lead]) - int(pos[p_lead])) / 1000,
        "TRAIT_AT_PROT_LEAD": float(t[p_lead]),
        "TRAIT_NEAR_PROT_LEAD": float(np.nanmax(t[near])),
        "PROT_AT_TRAIT_LEAD": float(p[t_lead]),
        "TRAIT_PEAKS": _clump(pos, t, 4.0),
        "PROT_PEAKS": _clump(pos, p, 5.0),
    }


def main() -> None:
    pl.Config.set_tbl_rows(200)
    pl.Config.set_tbl_width_chars(250)
    pl.Config.set_float_precision(2)
    rows = []
    for run, group in RUNS.items():
        out_dir = DEFAULT_RUNNER.meta_to_path(group.coloc_task.meta)
        raw = pl.read_parquet(out_dir / RAW_RESULTS_FILENAME).filter(
            pl.col(STATUS_COL) == STATUS_TESTED
        )
        for r in raw.iter_rows(named=True):
            stem = pair_file_stem(r[LOCUS_ID_COL], r[PROTEIN_ID_COL])
            window = pl.read_parquet(out_dir / WINDOWS_DIRNAME / f"{stem}.parquet")
            rows.append(
                {
                    "RUN": run,
                    LOCUS_ID_COL: r[LOCUS_ID_COL],
                    GENE_COL: r[GENE_COL],
                    PROTEIN_ID_COL: r[PROTEIN_ID_COL],
                    CIS_STATUS_COL: r[CIS_STATUS_COL],
                    PP_H3_COL: r[PP_H3_COL],
                    PP_H4_COL: r[PP_H4_COL],
                    **_pair_stats(window),
                }
            )
    df = pl.DataFrame(rows).unique([LOCUS_ID_COL, PROTEIN_ID_COL], keep="last")
    print(f"tested pairs (deduplicated across runs): {df.height}")
    print("\n== H3-dominant pairs, by trait signal under the protein lead (+/-25 kb)")
    print(
        df.filter(pl.col(PP_H3_COL) > 0.5).sort(
            "TRAIT_NEAR_PROT_LEAD", descending=True
        )
    )
    print(
        "\n== asymmetric pairs: protein strong at the trait lead, or trait strong near"
        " the protein lead"
    )
    print(
        df.filter(
            pl.col(PP_H3_COL) > 0.5,
            (pl.col("PROT_AT_TRAIT_LEAD") >= 4)
            | (pl.col("TRAIT_NEAR_PROT_LEAD") >= 4),
        )
        .select(
            "RUN",
            LOCUS_ID_COL,
            GENE_COL,
            CIS_STATUS_COL,
            PP_H3_COL,
            PP_H4_COL,
            "TRAIT_MAX",
            "PROT_MAX",
            "LEAD_DIST_KB",
            "TRAIT_AT_PROT_LEAD",
            "TRAIT_NEAR_PROT_LEAD",
            "PROT_AT_TRAIT_LEAD",
        )
        .sort("PROT_AT_TRAIT_LEAD", descending=True)
    )
    print("\n== loci with >1 distance-clumped trait peak")
    print(
        df.filter(pl.col("TRAIT_PEAKS") > 1)
        .group_by(LOCUS_ID_COL)
        .agg(pl.col("TRAIT_PEAKS").first(), pl.col("TRAIT_MAX").first(), pl.len())
        .sort("TRAIT_MAX", descending=True)
    )


if __name__ == "__main__":
    main()
