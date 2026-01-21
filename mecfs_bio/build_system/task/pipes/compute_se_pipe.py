"""
Compute standard error in GWAS summary statistics.
"""

import narwhals
import numpy as np
import pandas as pd
import scipy.stats

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_P_COL,
    GWASLAB_SE_COL,
)


class ComputeSEPipe(DataProcessingPipe):
    """
    Pipe to compute standard error in GWAS summary statistics using p-value and beta
    """

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        cols = x.collect_schema().keys()
        assert GWASLAB_SE_COL not in cols
        assert GWASLAB_P_COL in cols
        assert GWASLAB_BETA_COL in cols
        collected: pd.DataFrame = x.collect().to_pandas()
        pvals: np.ndarray = collected[GWASLAB_P_COL].to_numpy()
        z_score: np.ndarray = abs(scipy.stats.norm.ppf(1 - pvals / 2))
        z_score_2: np.ndarray = abs(scipy.stats.norm.ppf(pvals / 2))
        min_z_score = np.minimum(
            z_score, z_score_2
        )  # this is to resolve numerical issues when (1-pval) is numerically equal to 1
        se = abs(collected[GWASLAB_BETA_COL] / min_z_score)
        collected[GWASLAB_SE_COL] = se
        return narwhals.from_native(collected).lazy()


class ComputeSEIfNeededPipe(DataProcessingPipe):
    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        cols = x.collect_schema().keys()
        if GWASLAB_SE_COL not in cols:
            return ComputeSEPipe().process(x)
        return x
