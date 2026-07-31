"""
Pin down the mechanism: which per-SNP scaling of R's effect makes LDSC see the
beta-level orthogonality (gcov(R,pain)=0)?

Established so far:
  - Orthogonality holds at the genetic-covariance (beta) level: S12 - c*S22 = 0.
  - Feeding R's *kernel* output (Z_R = est/se_c, N=N_eff) into the SAME LDSC port
    gives rg(R,pain)=0.219 and an impossible h2(R)=3.07.
  - mean(N_eff)~2805 (input N's are 275k/388k): se_c is dominated by the global
    delta-method loading-uncertainty (V_LD) term, and corr(se_c,|est|)=+0.50.

A real single-trait GWAS has Z_j = beta_std,j * sqrt(N * varSNP_j) with N
CONSTANT. LDSC's recovery of genetic covariance assumes exactly that. The kernel
instead emits Z_R = est/se_c with a per-SNP, signal-correlated se_c -> the
constant-N assumption is violated.

TEST: rebuild R's sumstats with several Z_R definitions, all from the SAME betas
(est), and see which makes rg(R,pain) vanish:
  (a) kernel:    Z_R = Z_Estimate              (= est/se_c)   [pipeline baseline]
  (b) gwas-like: Z_R = est * sqrt(varSNP)       (constant-N GWAS scaling)
  (c) flat-se:   Z_R = est                      (se_c constant)
  (d) gwas-like but using the kernel N_eff column as the per-SNP N in LDSC, with
      Z_R = est*sqrt(varSNP) -- isolates whether N_eff (not Z) drives h2.

Whichever Z_R gives rg(R,pain)~0 identifies the scaling the pipeline should have
used, and hence the nature of the discrepancy.
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
LD_SRC = Path(
    "assets/base_asset_store/reference_data/linkage_disequilibrium_scores/"
    "thousand_genomes_phase_3_v1/extracted/"
    "thousand_genomes_phase_3_v1_eur_ld_scores_extracted"
)


def _strip(src: Path, dst: Path) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    for f in src.iterdir():
        if f.name.startswith("LDscore."):
            (dst / f.name[len("LDscore.") :]).symlink_to(f.resolve())


def _write(df: pl.DataFrame, z: np.ndarray, n: np.ndarray, out: Path) -> None:
    pl.DataFrame(
        {
            "SNP": df["SNP"],
            "N": n,
            "Z": z,
            "A1": df["A1"],
            "A2": df["A2"],
        }
    ).write_csv(out, separator="\t")


def _rg(tmpd: Path, ld_dir: Path, r_path: Path, label: str) -> None:
    res = run_ldsc(
        munged_paths=[MUNGED_D, MUNGED_PAIN, r_path],
        ld_dir=ld_dir,
        sample_prev=[0.0566, None, None],
        population_prev=[0.006, None, None],
    )
    ss = res.S_Stand
    print(f"\n[{label}]")
    print(f"  rg(R,pain) = {ss[1, 2]:+.4f}   rg(R,D) = {ss[0, 2]:+.4f}   "
          f"h2(R) = {res.S[2, 2]:.4f}   gcov(R,pain) = {res.S[1, 2]:+.5f}")


def main() -> None:
    df = pl.read_parquet(REMAINDER).filter(
        ~pl.col("fail")
        & pl.col("Z_Estimate").is_finite()
        & pl.col("N_eff").is_finite()
        & pl.col("est").is_finite()
        & pl.col("se_c").is_finite()
    )
    maf = df["MAF"].to_numpy()
    varsnp = 2 * maf * (1 - maf)
    est = df["est"].to_numpy()
    se_c = df["se_c"].to_numpy()
    z_kernel = df["Z_Estimate"].to_numpy()
    n_eff = df["N_eff"].to_numpy()

    # --- dissect se_c: how does it scale with varSNP? ---
    print("===== se_c dissection =====")
    print(f"se_c:  min={se_c.min():.4g} median={np.median(se_c):.4g} max={se_c.max():.4g}")
    # If se_c ~ 1/sqrt(N_eff*varSNP) (the N_eff definition), then se_c*sqrt(varSNP)
    # would be ~1/sqrt(N_eff). Check how flat se_c is vs varSNP scaling.
    print(f"corr(se_c, 1/sqrt(varSNP))   = {np.corrcoef(se_c, 1/np.sqrt(varsnp))[0,1]:+.3f}")
    print(f"corr(se_c, |est|)            = {np.corrcoef(se_c, np.abs(est))[0,1]:+.3f}")
    print(f"corr(se_c*sqrt(varSNP), |est|) = {np.corrcoef(se_c*np.sqrt(varsnp), np.abs(est))[0,1]:+.3f}")
    # check the identity N_eff == 1/(varSNP*se_c^2)
    n_check = 1.0 / (varsnp * se_c**2)
    print(f"max rel err of N_eff == 1/(varSNP*se_c^2): "
          f"{np.max(np.abs(n_check - n_eff) / n_eff):.2e}")

    with tempfile.TemporaryDirectory() as tmp:
        tmpd = Path(tmp)
        ld_dir = tmpd / "ld"
        _strip(LD_SRC, ld_dir)

        # (a) kernel baseline
        pa = tmpd / "a.sumstats"
        _write(df, z_kernel, n_eff, pa)
        _rg(tmpd, ld_dir, pa, "a) kernel: Z=Z_Estimate=est/se_c, N=N_eff")

        # (b) gwas-like scaling, constant N
        zb = est * np.sqrt(varsnp)
        zb = zb / np.std(zb)  # scale only; rg invariant
        pb = tmpd / "b.sumstats"
        _write(df, zb, np.full(len(zb), 10000.0), pb)
        _rg(tmpd, ld_dir, pb, "b) gwas-like: Z=est*sqrt(varSNP), N=const")

        # (c) flat se: Z proportional to beta
        zc = est / np.std(est)
        pc = tmpd / "c.sumstats"
        _write(df, zc, np.full(len(zc), 10000.0), pc)
        _rg(tmpd, ld_dir, pc, "c) flat-se: Z=est, N=const")

        # (d) gwas-like Z but with N_eff as the per-SNP N
        pd_ = tmpd / "d.sumstats"
        _write(df, zb, n_eff, pd_)
        _rg(tmpd, ld_dir, pd_, "d) gwas-like Z=est*sqrt(varSNP), N=N_eff")


if __name__ == "__main__":
    main()
    print("\nDONE", file=sys.stderr)
