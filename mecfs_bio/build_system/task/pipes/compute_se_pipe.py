"""
Compute standard error in GWAS summary statistics.
"""

import narwhals
import scipy.stats

from mecfs_bio.build_system.task.gwaslab.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_P_COL,
    GWASLAB_SE_COL,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


class ComputeSEPipe(DataProcessingPipe):
    """
    Pipe to compute standard error in GWAS summary statistics using p-value and beta
    """

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        cols = x.collect_schema().keys()
        assert GWASLAB_SE_COL not in cols
        assert GWASLAB_P_COL in cols
        assert GWASLAB_BETA_COL in cols
        collected = x.collect().to_pandas()
        z_score = scipy.stats.norm.ppf(1 - collected[GWASLAB_P_COL] / 2)
        se = abs(collected[GWASLAB_BETA_COL] / z_score)
        collected[GWASLAB_SE_COL] = se
        return narwhals.from_native(collected).lazy()
