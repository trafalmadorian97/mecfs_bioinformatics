import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_MLOG10P_COL, GWASLAB_P_COL


@frozen
class ComputeMlog10pIfNeededPipe(DataProcessingPipe):
    min_p_value: float = 1e-250

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        schema = x.collect_schema()
        if GWASLAB_MLOG10P_COL in schema:
            return x
        assert GWASLAB_P_COL in schema
        result = x.with_columns(
            (
                -1
                * (
                    narwhals.max_horizontal(
                        narwhals.col(GWASLAB_P_COL),
                        narwhals.lit(self.min_p_value),  # avoid taking the log of 0
                    ).log(
                        base=10,
                    )
                )
            ).alias(GWASLAB_MLOG10P_COL)
        )
        return result
