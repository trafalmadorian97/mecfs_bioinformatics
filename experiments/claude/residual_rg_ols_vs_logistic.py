"""
Final attribution test: is the residual rg(R,pain) caused specifically by the
binary/logistic composite trait (DecodeME)?

Mechanism claim: beta_R = (beta_D - c*beta_pain)/a_R is genetically orthogonal to
pain only if BOTH traits map to their LDSC Z the same way, Z = beta*sqrt(N*varSNP)
(the OLS map). Pain (OLS) does; DecodeME (logistic) does NOT -- its standardized
beta is logOR/sqrt(logOR^2*varSNP + pi^2/3), a different per-SNP map. So the
constant-c subtraction can't orthogonalize the logistic trait, but COULD
orthogonalize an OLS one.

TEST: reconstruct an OLS-style DecodeME beta from the munged Z (beta_D_ols =
Z_D/sqrt(N_D*varSNP)) and pain beta the same way, build R_ols(c) for a grid of c,
and show rg(R_ols, pain) passes cleanly through 0 -- i.e. an OLS-OLS subtraction
IS orthogonalizable by a constant c. Contrast with the real logistic composite,
which is stuck at rg=0.22 (established in earlier scripts).
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
MUNGED_D = ASSET / "munged" / "DECODE_ME.sumstats"
MUNGED_PAIN = ASSET / "munged" / "Multsite_Pain.sumstats"
REMAINDER = ASSET / "gwas_results" / "remainder_factor.parquet"
LD_SRC = Path(
    "assets/base_asset_store/reference_data/linkage_disequilibrium_scores/"
    "thousand_genomes_phase_3_v1/extracted/"
    "thousand_genomes_phase_3_v1_eur_ld_scores_extracted"
)
N_D, N_PAIN = 275488.0, 387649.0
C_KERNEL = 0.48399


def _strip(src: Path, dst: Path) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    for f in src.iterdir():
        if f.name.startswith("LDscore."):
            (dst / f.name[len("LDscore.") :]).symlink_to(f.resolve())


def main() -> None:
    # MAF (1000G) per SNP from the parquet.
    maf_df = pl.read_parquet(REMAINDER).select("SNP", "MAF", "A1", "A2").unique("SNP")
    d = pl.read_csv(MUNGED_D, separator="\t")[["SNP", "Z", "A1", "A2"]].rename(
        {"Z": "Zd", "A1": "A1d", "A2": "A2d"}
    )
    p = pl.read_csv(MUNGED_PAIN, separator="\t")[["SNP", "Z", "A1", "A2"]].rename(
        {"Z": "Zp", "A1": "A1p", "A2": "A2p"}
    )
    m = maf_df.join(d, on="SNP").join(p, on="SNP")
    # align everything to the parquet A1.
    a1 = m["A1"].to_numpy()
    zd = np.where(a1 == m["A1d"].to_numpy(), m["Zd"].to_numpy(), -m["Zd"].to_numpy())
    zp = np.where(a1 == m["A1p"].to_numpy(), m["Zp"].to_numpy(), -m["Zp"].to_numpy())
    maf = m["MAF"].to_numpy()
    varsnp = 2 * maf * (1 - maf)
    keep = (maf > 0) & (maf < 1)
    zd, zp, varsnp = zd[keep], zp[keep], varsnp[keep]
    snp = m["SNP"].to_numpy()[keep]
    A1 = m["A1"].to_numpy()[keep]
    A2 = m["A2"].to_numpy()[keep]
    n = len(snp)
    print(f"matched SNPs: {n}")

    # OLS-style standardized betas (same per-SNP map for BOTH traits).
    beta_d_ols = zd / np.sqrt(N_D * varsnp)
    beta_p = zp / np.sqrt(N_PAIN * varsnp)

    with tempfile.TemporaryDirectory() as tmp:
        tmpd = Path(tmp)
        ld_dir = tmpd / "ld"
        _strip(LD_SRC, ld_dir)

        def rg_for_c(c: float) -> tuple[float, float]:
            beta_r = beta_d_ols - c * beta_p  # a_R scaling irrelevant to rg
            zr = beta_r * np.sqrt(varsnp)
            zr = zr / np.std(zr)
            pr = tmpd / "r.sumstats"
            pl.DataFrame(
                {"SNP": snp, "N": 10000.0, "Z": zr, "A1": A1, "A2": A2}
            ).write_csv(pr, separator="\t")
            res = run_ldsc(
                munged_paths=[MUNGED_D, MUNGED_PAIN, pr],
                ld_dir=ld_dir,
                sample_prev=[0.0566, None, None],
                population_prev=[0.006, None, None],
            )
            ss = res.S_Stand
            return (ss[1, 2] if ss is not None else float("nan"), res.S[1, 2])

        print("\n===== OLS-OLS reconstruction: rg(R_ols, pain) vs subtraction coef c =====")
        print(f"{'c':>8} | {'rg(R_ols,pain)':>15} | {'gcov(R_ols,pain)':>16}")
        grid = [0.0, 0.2, 0.4, C_KERNEL, 0.55, 0.6, 0.7, 0.9]
        results = []
        for c in grid:
            rg, gcov = rg_for_c(c)
            results.append((c, rg))
            tag = "  <- kernel c" if abs(c - C_KERNEL) < 1e-3 else ""
            print(f"{c:8.4f} | {rg:15.4f} | {gcov:16.5f}{tag}")

        # crude root-find for the c giving rg=0 by linear interpolation.
        results.sort()
        cs = np.array([r[0] for r in results])
        rgs = np.array([r[1] for r in results])
        sign = np.sign(rgs)
        cross = np.where(np.diff(sign) != 0)[0]
        if len(cross):
            i = cross[0]
            c0 = cs[i] - rgs[i] * (cs[i + 1] - cs[i]) / (rgs[i + 1] - rgs[i])
            print(f"\n--> OLS-OLS subtraction is orthogonalizable: rg(R,pain)=0 at c ~ {c0:.4f}")
        print(f"    (the real LOGISTIC composite stays at rg=0.22 with kernel c={C_KERNEL};")
        print("     a constant c cannot zero it because Z_D's map to beta is non-OLS)")


if __name__ == "__main__":
    main()
    print("\nDONE", file=sys.stderr)
