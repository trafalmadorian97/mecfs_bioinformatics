"""
Probe the Kerrebijn fibromyalgia sumstats for the numerical pathology that makes
cross-trait LDSC fail with:

    ** On entry to DLASCL parameter number 4 had an illegal value
    numpy.linalg.LinAlgError: SVD did not converge in Linear Least Squares

That LAPACK message is emitted when the design matrix handed to lstsq contains
NaN/Inf, so this script looks for non-finite (or degenerate) values in the
columns that feed the regression: BETA, SE, N and the derived Z = BETA/SE.

It applies the same filtering the CT-LDSC task applies (indels, palindromes,
HLA, hapmap3, duplicate rsIDs) so the numbers describe the data that actually
reaches LDSC.

Run:
    pixi r python experiments/claude/fibro_ct_ldsc_nan_probe.py \
        2>&1 | tee experiments/claude/logs/fibro_ct_ldsc_nan_probe.log
"""

from pathlib import Path

import numpy as np
import pandas as pd

from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    FilterSettings,
    filter_sumstats,
)

FIBRO_PICKLE = Path(
    "assets/base_asset_store/gwas/fibromyalgia/kerrebijn_et_al/gwaslab_sumstats/kerrebijin_fibro_sumstats_37.pickle"
)
SCZ_PICKLE = Path(
    "assets/base_asset_store/gwas/schizophrenia/pgc_2022/gwaslab_sumstats/pgc_2022_sch_sumstats_37.pickle"
)


def describe(df: pd.DataFrame, label: str) -> None:
    print(f"\n===== {label} =====")
    print(f"rows: {len(df)}")
    print(f"columns: {list(df.columns)}")

    for col in ["BETA", "SE", "N", "EAF", "P"]:
        if col not in df.columns:
            print(f"  {col}: ABSENT")
            continue
        s = pd.to_numeric(df[col], errors="coerce")
        n_null = int(df[col].isna().sum())
        n_nan = int(np.isnan(s.to_numpy(dtype="float64")).sum())
        n_inf = int(np.isinf(s.to_numpy(dtype="float64")).sum())
        n_zero = int((s == 0).sum())
        n_neg = int((s < 0).sum())
        print(
            f"  {col}: dtype={df[col].dtype} null={n_null} nan={n_nan} inf={n_inf} "
            f"zero={n_zero} negative={n_neg} min={s.min()} max={s.max()}"
        )

    if "BETA" in df.columns and "SE" in df.columns:
        beta = pd.to_numeric(df["BETA"], errors="coerce").to_numpy(dtype="float64")
        se = pd.to_numeric(df["SE"], errors="coerce").to_numpy(dtype="float64")
        with np.errstate(divide="ignore", invalid="ignore"):
            z = beta / se
        n_nan = int(np.isnan(z).sum())
        n_inf = int(np.isinf(z).sum())
        print(
            f"  Z=BETA/SE: nan={n_nan} inf={n_inf} "
            f"finite_min={np.nanmin(z[np.isfinite(z)]) if np.isfinite(z).any() else 'NA'} "
            f"finite_max={np.nanmax(z[np.isfinite(z)]) if np.isfinite(z).any() else 'NA'}"
        )
        bad = ~np.isfinite(z)
        if bad.any():
            print(f"\n  --- first 10 rows with non-finite Z ({int(bad.sum())} total) ---")
            cols = [
                c
                for c in ["SNP", "rsID", "CHR", "POS", "EA", "NEA", "BETA", "SE", "P", "N", "EAF"]
                if c in df.columns
            ]
            print(df.loc[bad, cols].head(10).to_string())

        # chi2 = z^2 is what LDSC regresses; report the extremes
        chi2 = z[np.isfinite(z)] ** 2
        if chi2.size:
            print(
                f"  chi2: mean={chi2.mean():.4f} max={chi2.max():.4f} "
                f"n_over_80={int((chi2 > 80).sum())}"
            )

    if "N" in df.columns:
        n = pd.to_numeric(df["N"], errors="coerce")
        print(f"  N quantiles: {n.quantile([0, 0.001, 0.5, 0.999, 1.0]).to_dict()}")


def main() -> None:
    import gwaslab

    for label, path in [("FIBROMYALGIA", FIBRO_PICKLE), ("SCHIZOPHRENIA", SCZ_PICKLE)]:
        sumstats = gwaslab.load_pickle(str(path))
        describe(sumstats.data, f"{label} raw (as stored)")
        sumstats.infer_build()
        filter_sumstats(sumstats, FilterSettings(), build="19")
        describe(sumstats.data, f"{label} after CT-LDSC filtering")


if __name__ == "__main__":
    main()
