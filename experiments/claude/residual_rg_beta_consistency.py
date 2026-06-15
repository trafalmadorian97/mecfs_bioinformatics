"""
Decisive self-consistency test.

Established: rg(R,pain)=0.22 is intrinsic to R's betas (not se_c/N_eff/engine).
Since genetic covariance is LINEAR in per-SNP betas, and the kernel builds
    beta_R = (beta_D - c*beta_pain)/a_R,   c = S12/S22
orthogonality gcov(R,pain)=0 holds IFF, *measured on the betas the kernel
combined*, gcov(beta_D, beta_pain)/h2(beta_pain) == c. The catch: c comes from
ldsc() (munge pipeline, HapMap3), but the betas come from sumstats() (1000G).

This script recovers the EXACT sumstats betas the kernel used, without re-fetching
raw data:
  - beta_pain (OLS) is exactly reconstructable: beta_pain = Z_pain/sqrt(N_pain*varSNP)
  - beta_D is then recovered algebraically: beta_D = a_R*est + c*beta_pain
Then it runs LDSC on these betas (consistent GWAS scaling Z = beta*sqrt(varSNP*N))
and asks: do the betas REPRODUCE the kernel's S = [[.0826,.034],[.034,.0703]]?
If h2(beta_D) != 0.0826 or gcov(beta_D,beta_pain) != 0.034, the sumstats betas are
inconsistent with the ldsc S that set c -- that inconsistency is the root cause.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
import polars as pl

from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_ldsc import run_ldsc

ASSET = Path(
    "assets/base_asset_store/gwas/multi_trait/genomic_sem/analysis/"
    "decode_me_minus_pain_subtraction_full_python"
)
MUNGED_PAIN = ASSET / "munged" / "Multsite_Pain.sumstats"
REMAINDER = ASSET / "gwas_results" / "remainder_factor.parquet"
LD_SRC = Path(
    "assets/base_asset_store/reference_data/linkage_disequilibrium_scores/"
    "thousand_genomes_phase_3_v1/extracted/"
    "thousand_genomes_phase_3_v1_eur_ld_scores_extracted"
)

# kernel S (ldsc_S.csv) and derived loadings.
S11, S12, S22 = 0.08258, 0.034, 0.07025
C = S12 / S22
A_R = float(np.sqrt(S11 - S12**2 / S22))
N_D, N_PAIN = 275488.0, 387649.0


def _strip(src: Path, dst: Path) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    for f in src.iterdir():
        if f.name.startswith("LDscore."):
            (dst / f.name[len("LDscore.") :]).symlink_to(f.resolve())


def _write(snp, z, n, a1, a2, out: Path) -> None:
    pl.DataFrame({"SNP": snp, "N": n, "Z": z, "A1": a1, "A2": a2}).write_csv(
        out, separator="\t"
    )


def main() -> None:
    print(f"c = {C:.5f}   a_R = {A_R:.5f}")
    r = pl.read_parquet(REMAINDER).filter(
        ~pl.col("fail") & pl.col("est").is_finite() & pl.col("MAF").is_finite()
    ).select("SNP", "MAF", "A1", "A2", "est")
    pain = pl.read_csv(MUNGED_PAIN, separator="\t")[["SNP", "N", "Z", "A1", "A2"]]
    pain = pain.rename({"A1": "A1_p", "A2": "A2_p", "Z": "Z_p", "N": "N_p"})

    m = r.join(pain, on="SNP", how="inner")
    # align pain to the parquet (sumstats) A1.
    a1 = m["A1"].to_numpy()
    a1p = m["A1_p"].to_numpy()
    a2p = m["A2_p"].to_numpy()
    keep = (a1 == a1p) | (a1 == a2p)
    m = m.filter(pl.Series(keep))
    a1 = m["A1"].to_numpy(); a1p = m["A1_p"].to_numpy()
    z_p = np.where(a1 == a1p, m["Z_p"].to_numpy(), -m["Z_p"].to_numpy())

    maf = m["MAF"].to_numpy()
    varsnp = 2 * maf * (1 - maf)
    est = m["est"].to_numpy()

    # exact sumstats betas:
    beta_pain = z_p / np.sqrt(N_PAIN * varsnp)   # OLS reconstruction
    beta_D = A_R * est + C * beta_pain           # algebraic recovery of kernel input
    beta_R = est

    snp = m["SNP"]; A1 = m["A1"]; A2 = m["A2"]
    n_snp = m.height
    print(f"matched SNPs: {n_snp}")

    with tempfile.TemporaryDirectory() as tmp:
        tmpd = Path(tmp)
        ld_dir = tmpd / "ld"
        _strip(LD_SRC, ld_dir)

        # GWAS-consistent Z for each beta: Z = beta * sqrt(varSNP * N).
        pD = tmpd / "D.sumstats"
        _write(snp, beta_D * np.sqrt(varsnp * N_D), np.full(n_snp, N_D), A1, A2, pD)
        pP = tmpd / "P.sumstats"
        _write(snp, beta_pain * np.sqrt(varsnp * N_PAIN), np.full(n_snp, N_PAIN), A1, A2, pP)
        pR = tmpd / "R.sumstats"
        # R uses a nominal constant N (rg is N-invariant).
        _write(snp, beta_R * np.sqrt(varsnp), np.full(n_snp, 10000.0), A1, A2, pR)

        res = run_ldsc(
            munged_paths=[pD, pP, pR],
            ld_dir=ld_dir,
            sample_prev=[0.0566, None, None],
            population_prev=[0.006, None, None],
        )
        np.set_printoptions(precision=5, suppress=True)
        S = res.S
        print("\n===== LDSC on the RECONSTRUCTED sumstats betas =====")
        print("S_recon:\n", S)
        print("\nConsistency vs kernel S (the matrix that set c):")
        print(f"  h2(beta_D)        = {S[0,0]:.5f}   kernel S11 = {S11}")
        print(f"  h2(beta_pain)     = {S[1,1]:.5f}   kernel S22 = {S22}")
        print(f"  gcov(beta_D,pain) = {S[0,1]:.5f}   kernel S12 = {S12}")
        c_betas = S[0, 1] / S[1, 1]
        print(f"\n  c implied by sumstats betas = gcov(D,pain)/h2(pain) = {c_betas:.5f}")
        print(f"  c used by the kernel (S12/S22 from ldsc)            = {C:.5f}")
        print(f"  --> mismatch in c = {c_betas - C:+.5f}")
        if res.S_Stand is not None:
            ss = res.S_Stand
            print(f"\n  rg(R, pain) = {ss[1,2]:+.4f}   rg(R, D) = {ss[0,2]:+.4f}")
        print(f"  gcov(R,pain) = {S[1,2]:+.5f}  (should be 0 if betas satisfied c)")
        # predicted residual from the c-mismatch:
        print(f"\n  predicted gcov(R,pain) = (1/a_R)*h2(pain)*(c_betas - c) = "
              f"{(1/A_R)*S[1,1]*(c_betas - C):+.5f}")


if __name__ == "__main__":
    main()
    print("\nDONE", file=sys.stderr)
