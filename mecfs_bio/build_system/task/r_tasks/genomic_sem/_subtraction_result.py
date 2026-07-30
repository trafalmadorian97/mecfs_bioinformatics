"""
rpy2-free helpers for packaging GWAS-by-subtraction kernel output into the
per-factor result tables written to disk. Shared by the R-backed and the
full-Python subtraction tasks so both emit an identical schema.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from attrs import frozen

from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_config import (
    LDSC_BP_COL,
    LDSC_CHR_COL,
    MUNGE_A1_COL,
    MUNGE_A2_COL,
    MUNGE_MAF_COL,
    MUNGE_SNP_COL,
    SUBTRACTION_EST_COL,
    SUBTRACTION_FAIL_COL,
    SUBTRACTION_LHS_COL,
    SUBTRACTION_N_EFF_COL,
    SUBTRACTION_OP_COL,
    SUBTRACTION_P_COL,
    SUBTRACTION_RHS_COL,
    SUBTRACTION_SE_COL,
    SUBTRACTION_Z_COL,
)

# Model-term values placed in the lhs/op/rhs columns (the SEM "<factor> ~ SNP").
_MODEL_OP = "~"
_MODEL_RHS = "SNP"


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
            MUNGE_SNP_COL: snps_df[MUNGE_SNP_COL].to_numpy(),
            LDSC_CHR_COL: snps_df.get(LDSC_CHR_COL, pd.Series([pd.NA] * n)).to_numpy(),
            LDSC_BP_COL: snps_df.get(LDSC_BP_COL, pd.Series([pd.NA] * n)).to_numpy(),
            MUNGE_MAF_COL: snps_df[MUNGE_MAF_COL].to_numpy(),
            MUNGE_A1_COL: snps_df.get(MUNGE_A1_COL, pd.Series([""] * n)).to_numpy(),
            MUNGE_A2_COL: snps_df.get(MUNGE_A2_COL, pd.Series([""] * n)).to_numpy(),
            SUBTRACTION_LHS_COL: [lhs] * n,
            SUBTRACTION_OP_COL: [_MODEL_OP] * n,
            SUBTRACTION_RHS_COL: [_MODEL_RHS] * n,
            SUBTRACTION_EST_COL: est,
            SUBTRACTION_SE_COL: se_c,
            SUBTRACTION_Z_COL: z,
            SUBTRACTION_P_COL: p,
            SUBTRACTION_N_EFF_COL: n_eff,
            SUBTRACTION_FAIL_COL: fail,
        }
    )
