"""
Shared input guard for the LD score regression tasks.

Both the single-trait heritability regression and the cross-trait genetic
correlation regression feed Z scores to the same underlying LDSC estimator, and
both abort on a non-finite Z, so the guard lives here rather than in either task.
"""

import numpy as np
import pandas as pd
import structlog

from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.ldsc_constants import LDSC_Z_COL

logger = structlog.get_logger()


def drop_variants_with_degenerate_z(data: pd.DataFrame) -> pd.DataFrame:
    """Drop variants whose LDSC Z score would be non-finite.

    gwaslab builds the LDSC Z score as BETA / SE (or uses an existing Z column). A
    variant with SE == 0 therefore yields an infinite Z (or NaN, when BETA is also 0,
    as happens when a source reports an odds ratio rounded to 1.00). deCODE summary
    statistics contain many such variants: odds ratios rounded to 1.00 give BETA == 0
    and SE == 0, and underflowed p-values give SE == 0 with a non-zero BETA. The
    harmonised GWAS Catalog release of the Kerrebijn fibromyalgia GWAS does the same,
    reporting BETA as the smallest normal double alongside SE == 0. A non-finite Z
    makes the IRWLS reweighting produce a non-finite design matrix, which aborts the
    underlying SVD.

    LDSC's own munge step drops these variants; neither estimate_h2_by_ldsc nor
    estimate_rg_by_ldsc does (unlike the stratified path, which caps chi-square
    unconditionally), so we drop them here rather than at the call site, since the
    requirement is intrinsic to these regressions. Variants are matched to gwaslab's
    Z-resolution order: prefer an existing Z column, otherwise derive the finiteness
    requirement from BETA and SE.
    """
    if LDSC_Z_COL in data.columns:
        keep = np.isfinite(data[LDSC_Z_COL])
    elif GWASLAB_BETA_COL in data.columns and GWASLAB_SE_COL in data.columns:
        keep = (
            np.isfinite(data[GWASLAB_BETA_COL])
            & np.isfinite(data[GWASLAB_SE_COL])
            & (data[GWASLAB_SE_COL] > 0)
        )
    else:
        return data
    n_dropped = int((~keep).sum())
    if n_dropped:
        logger.info(
            "Dropped variants with degenerate (non-finite) LDSC Z score",
            n_dropped=n_dropped,
            n_remaining=int(keep.sum()),
        )
    return data.loc[keep].copy()
