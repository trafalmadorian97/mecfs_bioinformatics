"""
Investigation: why does ct-ldsc say the GWAS-by-subtraction remainder factor R
still has rg = 0.2332 with multisite pain, when the subtraction model implies
rg(R, pain) = 0?

THEORY (worked out algebraically): the kernel makes R's per-SNP Z a fixed
linear combination of Z_DecodeME and Z_pain. Running that through LDSC, the
cross-trait genetic covariance with pain is

    gcov(R, pain) = (1/a_R) * ( gcov(D, pain) - c * h2(pain) ),  c = S12/S22

which is EXACTLY zero iff the LDSC engine doing the evaluation recovers the same
gcov(D, pain) and h2(pain) that the subtraction's own LDSC used to set
c = S12/S22 = 0.4836 (S12=0.0340, S22=0.0703 from ldsc_S.csv).

This script runs the *same* LDSC engine that built S (the GenomicSEM Python
port, run_ldsc) on three traits [DECODE_ME, Multisite_Pain, R], where R is
reconstructed from the kernel's own remainder_factor.parquet (native Z_Estimate,
N_eff). If rg(R, pain) ~ 0 here, the kernel is internally consistent and the
0.2332 is an engine/munge difference (gwaslab ct-ldsc vs the GenomicSEM port).
If rg(R, pain) ~ 0.23 even here, the leak is intrinsic to R's construction.

As a harness validation, we first reproduce the saved 2-trait S (D, pain) and
check S[0,1]=0.0340, S[1,1]=0.0703.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
import polars as pl

from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_ldsc import run_ldsc

ASSET = Path(
    "assets/base_asset_store/gwas/multi_trait/genomic_sem/analysis/"
    "decode_me_minus_pain_subtraction_full_python"
)
MUNGED_D = ASSET / "munged" / "DECODE_ME.sumstats"
MUNGED_PAIN = ASSET / "munged" / "Multsite_Pain.sumstats"
REMAINDER = ASSET / "gwas_results" / "remainder_factor.parquet"
COMMON = ASSET / "gwas_results" / "common_factor.parquet"

LD_SRC = Path(
    "assets/base_asset_store/reference_data/linkage_disequilibrium_scores/"
    "thousand_genomes_phase_3_v1/extracted/"
    "thousand_genomes_phase_3_v1_eur_ld_scores_extracted"
)

# DecodeME GWAS-1 prevalences (from prevalance_info.py) for liability scaling,
# so the 2-trait run reproduces the saved liability-scale S.
DECODE_SAMPLE_PREV = 0.0566
DECODE_POP_PREV = 0.006


def _strip_ldscore_prefix(src: Path, dst: Path) -> None:
    """Symlink LDscore.<chr>.* -> <chr>.* so run_ldsc's naming works."""
    dst.mkdir(parents=True, exist_ok=True)
    for f in src.iterdir():
        if f.name.startswith("LDscore."):
            (dst / f.name[len("LDscore.") :]).symlink_to(f.resolve())


def _write_factor_sumstats(parquet: Path, out: Path) -> int:
    """Build an LDSC .sumstats (SNP N Z A1 A2) from a kernel factor parquet."""
    df = pl.read_parquet(parquet)
    df = df.filter(
        ~pl.col("fail")
        & pl.col("Z_Estimate").is_finite()
        & pl.col("N_eff").is_finite()
    )
    out_df = df.select(
        pl.col("SNP"),
        pl.col("N_eff").alias("N"),
        pl.col("Z_Estimate").alias("Z"),
        pl.col("A1"),
        pl.col("A2"),
    )
    out_df.write_csv(out, separator="\t")
    return out_df.height


def _print_matrices(label: str, names: list[str], res) -> None:
    print(f"\n===== {label} =====")
    np.set_printoptions(precision=5, suppress=True)
    print("traits:", names)
    print("S (genetic covariance):\n", res.S)
    print("diag(S) = h2:", np.diag(res.S))
    if res.S_Stand is not None:
        print("S_Stand (genetic correlation):\n", res.S_Stand)
    print("I (intercepts):\n", res.I)


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpd = Path(tmp)
        ld_dir = tmpd / "ld"
        _strip_ldscore_prefix(LD_SRC, ld_dir)

        r_path = tmpd / "R.sumstats"
        n_r = _write_factor_sumstats(REMAINDER, r_path)
        print(f"R sumstats rows (post-filter): {n_r}")

        # ---- Harness validation: reproduce the saved 2-trait S ----
        res2 = run_ldsc(
            munged_paths=[MUNGED_D, MUNGED_PAIN],
            ld_dir=ld_dir,
            sample_prev=[DECODE_SAMPLE_PREV, None],
            population_prev=[DECODE_POP_PREV, None],
        )
        _print_matrices(
            "VALIDATION: 2-trait [DECODE_ME, pain] (should match ldsc_S.csv)",
            ["DECODE_ME", "pain"],
            res2,
        )
        print("Expected from ldsc_S.csv: S[0,0]=0.08258 S[0,1]=0.03400 S[1,1]=0.07025")
        c_internal = res2.S[0, 1] / res2.S[1, 1]
        print(f"c = S12/S22 (this run) = {c_internal:.5f}  (kernel used 0.4836)")

        # ---- Decisive run: 3-trait [DECODE_ME, pain, R] ----
        # R is treated as quantitative (no prevalence); rg is liability-invariant.
        res3 = run_ldsc(
            munged_paths=[MUNGED_D, MUNGED_PAIN, r_path],
            ld_dir=ld_dir,
            sample_prev=[DECODE_SAMPLE_PREV, None, None],
            population_prev=[DECODE_POP_PREV, None, None],
        )
        _print_matrices(
            "DECISIVE: 3-trait [DECODE_ME, pain, R] via GenomicSEM-port LDSC",
            ["DECODE_ME", "pain", "R"],
            res3,
        )
        ss = res3.S_Stand
        print("\n---- genetic correlations (GenomicSEM-port LDSC) ----")
        print(f"rg(DECODE_ME, pain) = {ss[0, 1]:.4f}   (ct-ldsc external: 0.4535)")
        print(f"rg(R,         pain) = {ss[1, 2]:.4f}   (ct-ldsc external: 0.2332)")
        print(f"rg(R,    DECODE_ME) = {ss[0, 2]:.4f}   (ct-ldsc external: 0.9748)")

        # Predicted residual from the algebra, using THIS engine's gcov/h2:
        gcov_dp = res3.S[0, 1]
        h2_pain = res3.S[1, 1]
        c = res2.S[0, 1] / res2.S[1, 1]
        a_r = float(np.sqrt(res2.S[0, 0] - (res2.S[0, 1] ** 2) / res2.S[1, 1]))
        gcov_r_pain_pred = (1.0 / a_r) * (gcov_dp - c * h2_pain)
        print("\n---- algebraic decomposition (same-engine, should ~cancel) ----")
        print(f"gcov(D,pain)={gcov_dp:.5f}  c*h2(pain)={c * h2_pain:.5f}  "
              f"=> gcov(R,pain) predicted = {gcov_r_pain_pred:.6e}")
        print(f"gcov(R,pain) measured (S[1,2]) = {res3.S[1, 2]:.6e}")


if __name__ == "__main__":
    main()
    print("\nDONE", file=sys.stderr)
