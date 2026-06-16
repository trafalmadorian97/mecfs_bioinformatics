"""
Follow-up: pin down WHY rg(R, pain) != 0 even in the same LDSC engine.

Finding from residual_rg_internal_ldsc.py:
  - Beta-level identity says gcov(R,pain) = (1/a_R)(S12 - c*S22) = 0 exactly.
  - But LDSC MEASURES gcov(R,pain) = 0.102 and h2(R) = 3.07 (impossible: R is a
    unit-variance latent factor, h2 should be <= ~0.08-ish).

Hypothesis: the kernel's per-SNP Z_R = beta_R / se_c is NOT a constant-coefficient
linear combination of Z_D and Z_pain, because se_c (delta-method sandwich SE)
includes the V_LD loading-uncertainty term, which is signal-dependent and
correlated across traits. That makes N_eff = 1/(varSNP*se_c^2) heterogeneous and
correlated with effect sizes, violating LDSC's constant-N assumption.

TEST: construct a CLEAN Z_R as the exact constant-coefficient combination
    Z_R_clean ∝ Z_D / sqrt(N_D) - c * Z_pain / sqrt(N_pain)
(which the per-SNP-only sandwich SE would give), and run the same LDSC. If
rg(R_clean, pain) ~ 0 and h2(R_clean) is sane, the kernel's se_c / N_eff is the
culprit (not the betas, not the loadings, not the engine).

Also: directly measure whether the kernel's N_eff is heterogeneous and
correlated with the per-SNP signal.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_ldsc import run_ldsc

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

# Loadings from the saved S (ldsc_S.csv): S = [[0.08258,0.034],[0.034,0.07025]]
S11, S12, S22 = 0.08258, 0.034, 0.07025
C = S12 / S22  # = a_F/b = 0.4836
A_R = float(np.sqrt(S11 - S12**2 / S22))
N_D = 275488.0
N_PAIN = 387649.0


def _strip_ldscore_prefix(src: Path, dst: Path) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    for f in src.iterdir():
        if f.name.startswith("LDscore."):
            (dst / f.name[len("LDscore.") :]).symlink_to(f.resolve())


def _build_clean_R(out: Path) -> int:
    """Z_R_clean = Z_D/sqrt(N_D) - C*Z_pain/sqrt(N_pain), alleles aligned to D's A1."""
    d = pd.read_csv(MUNGED_D, sep="\t")[["SNP", "N", "Z", "A1", "A2"]]
    p = pd.read_csv(MUNGED_PAIN, sep="\t")[["SNP", "N", "Z", "A1", "A2"]]
    m = d.merge(p, on="SNP", suffixes=("_d", "_p"))
    # align pain Z to D's A1 (simple allele-pair agreement; both are munged to a
    # common reference so A1 either matches or is swapped).
    same = (m["A1_d"] == m["A1_p"]).to_numpy()
    swapped = (m["A1_d"] == m["A2_p"]).to_numpy()
    keep = same | swapped
    m = m.loc[keep].reset_index(drop=True)
    same = (m["A1_d"] == m["A1_p"]).to_numpy()
    z_p_aligned = np.where(same, m["Z_p"].to_numpy(), -m["Z_p"].to_numpy())
    z_r = m["Z_d"].to_numpy() / np.sqrt(N_D) - C * z_p_aligned / np.sqrt(N_PAIN)
    # scale to unit variance (constant scale; irrelevant to rg). Use a constant N.
    z_r = z_r / np.std(z_r)
    out_df = pd.DataFrame(
        {"SNP": m["SNP"], "N": 10000.0, "Z": z_r, "A1": m["A1_d"], "A2": m["A2_d"]}
    )
    out_df.to_csv(out, sep="\t", index=False)
    return len(out_df)


def _neff_diagnostics() -> None:
    print("\n===== N_eff / se_c heterogeneity diagnostics (kernel R parquet) =====")
    df = pl.read_parquet(REMAINDER).filter(
        ~pl.col("fail")
        & pl.col("N_eff").is_finite()
        & pl.col("Z_Estimate").is_finite()
    )
    neff = df["N_eff"].to_numpy()
    print(f"N_eff: min={neff.min():.0f}  median={np.median(neff):.0f}  "
          f"mean={neff.mean():.0f}  max={neff.max():.0f}  "
          f"CV={neff.std() / neff.mean():.3f}")
    # Is N_eff correlated with the per-SNP signal? In a normal GWAS N is ~constant
    # per SNP. If N_eff varies with MAF / effect size, LDSC's constant-N
    # assumption is violated.
    maf = df["MAF"].to_numpy()
    est = df["est"].to_numpy()
    se = df["se_c"].to_numpy()
    varsnp = 2 * maf * (1 - maf)
    print(f"corr(N_eff, MAF)            = {np.corrcoef(neff, maf)[0, 1]:+.3f}")
    print(f"corr(N_eff, varSNP)         = {np.corrcoef(neff, varsnp)[0, 1]:+.3f}")
    print(f"corr(N_eff, est^2)          = {np.corrcoef(neff, est**2)[0, 1]:+.3f}")
    print(f"corr(se_c, |est|)           = {np.corrcoef(se, np.abs(est))[0, 1]:+.3f}")
    print("(a normal GWAS would have ~constant N per SNP; strong corr with MAF or "
          "effect size signals heterogeneous, signal-dependent N_eff)")


def main() -> None:
    print(f"C = S12/S22 = {C:.5f}   a_R = {A_R:.5f}")
    with tempfile.TemporaryDirectory() as tmp:
        tmpd = Path(tmp)
        ld_dir = tmpd / "ld"
        _strip_ldscore_prefix(LD_SRC, ld_dir)

        rclean = tmpd / "R_clean.sumstats"
        n = _build_clean_R(rclean)
        print(f"R_clean SNPs: {n}")

        res = run_ldsc(
            munged_paths=[MUNGED_D, MUNGED_PAIN, rclean],
            ld_dir=ld_dir,
            sample_prev=[0.0566, None, None],
            population_prev=[0.006, None, None],
        )
        np.set_printoptions(precision=5, suppress=True)
        print("\n===== [D, pain, R_clean] (constant-coeff combination) =====")
        print("S:\n", res.S)
        print("diag(S)=h2:", np.diag(res.S))
        ss = res.S_Stand
        print("S_Stand:\n", ss)
        print(f"\nrg(D, pain)      = {ss[0, 1]:.4f}")
        print(f"rg(R_clean, pain)= {ss[1, 2]:.4f}   <-- expect ~0 if mechanism is se_c/N_eff")
        print(f"rg(R_clean, D)   = {ss[0, 2]:.4f}")
        print(f"h2(R_clean)      = {res.S[2, 2]:.4f}   <-- expect sane (<~0.1), not 3.07")

    _neff_diagnostics()


if __name__ == "__main__":
    main()
    print("\nDONE", file=sys.stderr)
