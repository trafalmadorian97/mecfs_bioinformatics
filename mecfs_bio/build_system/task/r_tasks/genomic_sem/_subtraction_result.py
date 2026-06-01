"""
rpy2-free helpers for packaging GWAS-by-subtraction kernel output into the
per-factor result tables written to disk. Shared by the R-backed and the
full-Python subtraction tasks so both emit an identical schema.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from attrs import frozen


@frozen
class SubtractionFrames:
    """The two per-factor result tables written by a subtraction task."""

    f_df: pd.DataFrame  # common factor F ~ SNP
    r_df: pd.DataFrame  # remainder factor R ~ SNP


def make_result_df(
    snps_df: pd.DataFrame,
    est: np.ndarray,
    se_c: np.ndarray,
    z: np.ndarray,
    p: np.ndarray,
    n_eff: np.ndarray,
    fail: np.ndarray,
    lhs: str,
) -> pd.DataFrame:
    n = len(snps_df)
    return pd.DataFrame(
        {
            "SNP": snps_df["SNP"].to_numpy(),
            "CHR": snps_df.get("CHR", pd.Series([pd.NA] * n)).to_numpy(),
            "BP": snps_df.get("BP", pd.Series([pd.NA] * n)).to_numpy(),
            "MAF": snps_df["MAF"].to_numpy(),
            "A1": snps_df.get("A1", pd.Series([""] * n)).to_numpy(),
            "A2": snps_df.get("A2", pd.Series([""] * n)).to_numpy(),
            "lhs": [lhs] * n,
            "op": ["~"] * n,
            "rhs": ["SNP"] * n,
            "est": est,
            "se_c": se_c,
            "Z_Estimate": z,
            "Pval_Estimate": p,
            "N_eff": n_eff,
            "fail": fail,
        }
    )
