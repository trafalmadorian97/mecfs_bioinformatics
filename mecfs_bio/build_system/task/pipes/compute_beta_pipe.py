import narwhals

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_ODDS_RATIO_COL,
)


class ComputeBetaPipe(DataProcessingPipe):
    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        schema = x.collect_schema()
        assert GWASLAB_BETA_COL not in schema.keys()
        x = x.with_columns(
            narwhals.col(GWASLAB_ODDS_RATIO_COL).log().alias(GWASLAB_BETA_COL)
        )
        return x
