"""
Compute standard error in GWAS summary statistics.
"""
from typing import Literal

import narwhals
import numpy as np
import pandas as pd
import polars
import scipy.stats
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_P_COL,
    GWASLAB_SE_COL,
)

Backend = Literal["polars"]

@frozen
class ComputeSEPipe(DataProcessingPipe):
    """
    Pipe to compute standard error in GWAS summary statistics using p-value and beta
    """
    min_p_value: float=1e-250
    max_p_value: float=0.99999999
    min_se_value: float=1e-20
    out_backend: None| Backend = None


    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        cols = x.collect_schema().keys()
        assert GWASLAB_SE_COL not in cols
        assert GWASLAB_P_COL in cols
        assert GWASLAB_BETA_COL in cols
        collected: pd.DataFrame = x.collect().to_pandas()
        pvals: np.ndarray = collected[GWASLAB_P_COL].to_numpy()
        pvals =np.maximum(pvals, self.min_p_value)
        pvals =np.minimum(pvals, self.max_p_value)
        z_score: np.ndarray = abs(scipy.stats.norm.ppf(1 - pvals / 2))
        z_score_2: np.ndarray = abs(scipy.stats.norm.ppf(pvals / 2))
        min_z_score = np.minimum(
            z_score, z_score_2
        )  # this is to resolve numerical issues when (1-pval) is numerically equal to 1
        se = abs(collected[GWASLAB_BETA_COL] / min_z_score)
        se = np.maximum(se, self.min_se_value)
        collected[GWASLAB_SE_COL] = se
        if self.out_backend is None:
            return narwhals.from_native(collected).lazy()
        if self.out_backend == "polars":
            return narwhals.from_native(polars.from_pandas(collected)).lazy()


class ComputeSEIfNeededPipe(DataProcessingPipe):
    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        cols = x.collect_schema().keys()
        if GWASLAB_SE_COL not in cols:
            return ComputeSEPipe().process(x)
        return x
