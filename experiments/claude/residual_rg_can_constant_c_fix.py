"""
Can a CONSTANT subtraction coefficient orthogonalize the REAL (logistic) R?

This decides whether the binary-composite scale mismatch is:
  (i) a global scalar  -> some c' != 0.484 fully zeroes rg(R,pain). Fix = "use the
      right c" (compute c on the same scale you subtract on).
  (ii) per-SNP structural (logistic den = sqrt(logOR^2*varSNP+pi^2/3) varies with
      effect size) -> rg bottoms out ABOVE 0; no constant c works, need a different
      formulation.

Method: recover the exact logistic beta_D the kernel used from its own output:
    beta_D = a_R*est + c_kernel*beta_pain      (since est=(beta_D - c*beta_pain)/a_R)
with beta_pain reconstructed exactly from the munged pain (OLS). Then sweep c' in
R(c') = beta_D - c'*beta_pain, scale Z=R*sqrt(varSNP) (proper GWAS map, benign
chisq), and run LDSC vs the real munged pain.
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
S11, S12, S22 = 0.08258, 0.034, 0.07025
C_KERNEL = S12 / S22
A_R = float(np.sqrt(S11 - S12**2 / S22))
N_PAIN = 387649.0


def _strip(src: Path, dst: Path) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    for f in src.iterdir():
        if f.name.startswith("LDscore."):
            (dst / f.name[len("LDscore.") :]).symlink_to(f.resolve())


def main() -> None:
    print(f"c_kernel = {C_KERNEL:.4f}   a_R = {A_R:.4f}")
    r = pl.read_parquet(REMAINDER).filter(
        ~pl.col("fail") & pl.col("est").is_finite()
    ).select("SNP", "MAF", "A1", "A2", "est")
    p = pl.read_csv(MUNGED_PAIN, separator="\t")[["SNP", "Z", "A1", "A2"]].rename(
        {"Z": "Zp", "A1": "A1p", "A2": "A2p"}
    )
    m = r.join(p, on="SNP", how="inner")
    a1 = m["A1"].to_numpy()
    keep = (a1 == m["A1p"].to_numpy()) | (a1 == m["A2p"].to_numpy())
    m = m.filter(pl.Series(keep))
    a1 = m["A1"].to_numpy()
    zp = np.where(a1 == m["A1p"].to_numpy(), m["Zp"].to_numpy(), -m["Zp"].to_numpy())
    maf = m["MAF"].to_numpy()
    varsnp = 2 * maf * (1 - maf)
    est = m["est"].to_numpy()
    snp, A1, A2 = m["SNP"].to_numpy(), m["A1"].to_numpy(), m["A2"].to_numpy()
    nrow = m.height

    beta_pain = zp / np.sqrt(N_PAIN * varsnp)        # OLS, exact
    beta_D = A_R * est + C_KERNEL * beta_pain        # exact recovery of kernel's logistic beta_D
    print(f"matched SNPs: {nrow}")

    with tempfile.TemporaryDirectory() as tmp:
        tmpd = Path(tmp)
        ld_dir = tmpd / "ld"
        _strip(LD_SRC, ld_dir)
        pain_path = MUNGED_PAIN

        def rg_for_c(cp: float) -> float:
            beta_r = beta_D - cp * beta_pain
            zr = beta_r * np.sqrt(varsnp)
            zr = zr / np.std(zr)
            pr = tmpd / "r.sumstats"
            pl.DataFrame(
                {"SNP": snp, "N": 10000.0, "Z": zr, "A1": A1, "A2": A2}
            ).write_csv(pr, separator="\t")
            res = run_ldsc(
                munged_paths=[pain_path, pr],
                ld_dir=ld_dir,
                sample_prev=[None, None],
                population_prev=[None, None],
            )
            return res.S_Stand[0, 1] if res.S_Stand is not None else float("nan")

        print("\n===== REAL logistic R: rg(R,pain) vs subtraction coef c' =====")
        print(f"{'c-prime':>8} | {'rg(R,pain)':>11}")
        grid = [0.0, 0.484, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5]
        pts = []
        for cp in grid:
            rg = rg_for_c(cp)
            pts.append((cp, rg))
            tag = "  <- kernel c" if abs(cp - C_KERNEL) < 1e-3 else ""
            print(f"{cp:8.3f} | {rg:11.4f}{tag}")

        pts.sort()
        cs = np.array([x[0] for x in pts]); rgs = np.array([x[1] for x in pts])
        cross = np.where(np.diff(np.sign(rgs)) != 0)[0]
        if len(cross):
            i = cross[0]
            c0 = cs[i] - rgs[i] * (cs[i + 1] - cs[i]) / (rgs[i + 1] - rgs[i])
            rg_min = rg_for_c(c0)
            print(f"\n--> rg(R,pain) crosses 0 near c' = {c0:.3f} (rg there = {rg_min:+.4f})")
            print(f"    kernel used c = {C_KERNEL:.3f}; rg at kernel c = {rg_for_c(C_KERNEL):+.4f}")
            print("    If a constant c' cleanly zeroes it -> global-scalar fix (recompute c "
                  "on the subtraction scale). Residual at the crossing = the per-SNP "
                  "(logistic-den) part that no constant c can remove.")
        else:
            print("\n--> rg never crosses 0 over the grid: per-SNP structural mismatch dominates.")


if __name__ == "__main__":
    main()
    print("\nDONE", file=sys.stderr)
