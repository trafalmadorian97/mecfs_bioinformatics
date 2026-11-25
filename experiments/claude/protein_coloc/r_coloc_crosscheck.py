"""
Cross-checks our Python coloc.abf (ColocAbf) against the reference R implementation (coloc
5.2.3) on the exploratory DecodeME x UKB-PPP coloc run.

For every tested (locus, protein) pair, the aligned window the coloc task saved (the exact
variants and betas ColocAbf used) is exported to CSV and run through R coloc::coloc.abf with the
same priors (P1, P2, P12 from raw_results) and the same trait models: the trait as case-control
(prior SD 0.2) and the protein as quantitative with sdY = 1 (prior SD 0.15). R runs in a
throwaway pixi exec environment, so coloc is not added to the project dependencies.

This verifies the Bayes factor and posterior arithmetic only, not the harmonization and
alignment that produced the windows.

Usage: pixi r python experiments/claude/protein_coloc/r_coloc_crosscheck.py
Exits non-zero when any posterior differs by 1e-6 or more.
"""

import sys
import tempfile
from pathlib import Path

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_ppp_coloc import (
    DECODE_ME_PPP_COLOC_EXPLORATORY_2,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import CaseControlTraitModel
from mecfs_bio.constants.protein_coloc_constants import (
    BETA_PROTEIN_COL,
    BETA_TRAIT_COL,
    GENE_COL,
    LOCUS_ID_COL,
    N_ALIGNED_COL,
    PP_H0_COL,
    PP_H1_COL,
    PP_H2_COL,
    PP_H3_COL,
    PP_H4_COL,
    PRIOR_P1_COL,
    PRIOR_P2_COL,
    PRIOR_P12_COL,
    PROTEIN_ID_COL,
    RAW_RESULTS_FILENAME,
    SE_PROTEIN_COL,
    SE_TRAIT_COL,
    STATUS_COL,
    STATUS_TESTED,
    WINDOWS_DIRNAME,
    pair_file_stem,
)
from mecfs_bio.util.subproc.run_command import execute_command

TOLERANCE = 1e-6
HIT_THRESHOLD = 0.8
R_SCRIPT = Path(__file__).with_name("r_coloc_crosscheck.R")
PP_COLS = [PP_H0_COL, PP_H1_COL, PP_H2_COL, PP_H3_COL, PP_H4_COL]


def export_windows(tested: pl.DataFrame, coloc_dir: Path, out_dir: Path, s: float) -> Path:
    manifest_rows = []
    for row in tested.iter_rows(named=True):
        stem = pair_file_stem(row[LOCUS_ID_COL], row[PROTEIN_ID_COL])
        window = pl.read_parquet(coloc_dir / WINDOWS_DIRNAME / f"{stem}.parquet")
        csv_path = out_dir / f"{stem}.csv"
        window.select(
            pl.format("{}:{}:{}:{}", "CHR", "POS", "EA", "NEA").alias("snp"),
            BETA_TRAIT_COL,
            SE_TRAIT_COL,
            BETA_PROTEIN_COL,
            SE_PROTEIN_COL,
        ).write_csv(csv_path)
        manifest_rows.append(
            {
                "pair": stem,
                "path": str(csv_path),
                "p1": row[PRIOR_P1_COL],
                "p2": row[PRIOR_P2_COL],
                "p12": row[PRIOR_P12_COL],
                "s": s,
            }
        )
    manifest = out_dir / "manifest.csv"
    pl.DataFrame(manifest_rows).write_csv(manifest)
    return manifest


def run_r(manifest: Path, out_csv: Path) -> pl.DataFrame:
    execute_command(
        [
            "pixi",
            "exec",
            "--spec",
            "r-base=4.4",
            "--spec",
            "r-coloc=5.2.3",
            "Rscript",
            str(R_SCRIPT),
            str(manifest),
            str(out_csv),
        ]
    )
    return pl.read_csv(out_csv)


def report(compared: pl.DataFrame, label: str) -> float:
    diffs = compared.select(
        [(pl.col(c) - pl.col(f"R_{c}")).abs().max().alias(c) for c in PP_COLS]
    )
    worst = max(diffs.row(0))
    print(f"--- {label}: {compared.height} pairs; max |dPP| per hypothesis:")
    print(diffs)
    print(f"overall max |dPP| = {worst:.3e}")
    return worst


def main() -> int:
    coloc_task = DECODE_ME_PPP_COLOC_EXPLORATORY_2.coloc_task
    trait_model = coloc_task.trait_model
    assert isinstance(trait_model, CaseControlTraitModel)
    coloc_dir = Path(DEFAULT_RUNNER.meta_to_path(coloc_task.meta))
    raw = pl.read_parquet(coloc_dir / RAW_RESULTS_FILENAME)
    tested = raw.filter(pl.col(STATUS_COL) == STATUS_TESTED)
    print(f"coloc dir: {coloc_dir}; tested pairs: {tested.height}")
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        manifest = export_windows(
            tested, coloc_dir, tmp_dir, s=trait_model.case_fraction
        )
        r_out = run_r(manifest, tmp_dir / "r_out.csv")
    keyed = tested.with_columns(
        pl.struct(LOCUS_ID_COL, PROTEIN_ID_COL)
        .map_elements(
            lambda r: pair_file_stem(r[LOCUS_ID_COL], r[PROTEIN_ID_COL]),
            return_dtype=pl.String,
        )
        .alias("pair")
    )
    compared = keyed.join(r_out, on="pair", how="inner")
    assert compared.height == tested.height, (
        f"R returned {r_out.height} pairs for {tested.height} tested"
    )
    n_mismatch = compared.filter(pl.col("nsnps") != pl.col(N_ALIGNED_COL)).height
    print(f"pairs where R nsnps != N_ALIGNED: {n_mismatch}")
    worst_all = report(compared, "all tested pairs")
    report(compared.filter(pl.col(PP_H4_COL) > HIT_THRESHOLD), "hits (PP.H4 > 0.8)")
    print("--- 5 pairs with the largest |dPP.H4|:")
    print(
        compared.with_columns(
            (pl.col(PP_H4_COL) - pl.col(f"R_{PP_H4_COL}")).abs().alias("dPP_H4")
        )
        .sort("dPP_H4", descending=True)
        .head(5)
        .select(LOCUS_ID_COL, GENE_COL, PP_H4_COL, f"R_{PP_H4_COL}", "dPP_H4")
    )
    print("--- hits, Python vs R PP.H4:")
    print(
        compared.filter(pl.col(PP_H4_COL) > HIT_THRESHOLD)
        .sort(PP_H4_COL, descending=True)
        .select(GENE_COL, LOCUS_ID_COL, PP_H4_COL, f"R_{PP_H4_COL}", "nsnps")
    )
    ok = worst_all < TOLERANCE and n_mismatch == 0
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
