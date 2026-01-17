import narwhals

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_MLOG10P_COL,
    GWASLAB_P_COL,
)


class ComputePPipe(DataProcessingPipe):
    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        schema = x.collect_schema()
        assert GWASLAB_P_COL not in schema
        assert GWASLAB_MLOG10P_COL in schema
        return x.with_columns(
            (10 ** (-1 * narwhals.col(GWASLAB_MLOG10P_COL))).alias(GWASLAB_P_COL)
        )


class ComputePIfNeededPipe(DataProcessingPipe):
    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        schema = x.collect_schema()
        if GWASLAB_P_COL in schema:
            return x
        return ComputePPipe().process(x)
