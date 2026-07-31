import narwhals
import polars as pl
import scipy.stats
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_P_COL,
    GWASLAB_SE_COL,
)


@frozen
class ComputePFromBetaSEPipeIfNeeded(DataProcessingPipe):
    p_col: str = GWASLAB_P_COL
    se_col: str = GWASLAB_SE_COL
    beta_col: str = GWASLAB_BETA_COL

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        schema = x.collect_schema()
        if self.p_col in schema:
            return x
        assert self.se_col in schema
        assert self.beta_col in schema
        collected = x.collect().to_pandas()
        z = collected[self.beta_col] / collected[self.se_col]
        collected[self.p_col] = 2 * scipy.stats.norm.sf(abs(z))
        collected_polars = pl.from_pandas(collected)
        return narwhals.from_native(collected_polars).lazy()
