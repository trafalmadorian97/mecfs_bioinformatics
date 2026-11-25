import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_95L_ODD_RATIO_COL,
    GWASLAB_95U_ODD_RATIO_COL,
    GWASLAB_SE_COL,
)


@frozen(slots=True)
class ComputeSEFromOddsCiPipe(DataProcessingPipe):
    """
    NOTE: GWASLAB SEEMS TO MESS UP LIFTOVER
    """

    min_se: float = 0

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        schema = x.collect_schema()
        assert GWASLAB_SE_COL not in schema
        assert GWASLAB_95L_ODD_RATIO_COL in schema
        assert GWASLAB_95U_ODD_RATIO_COL in schema
        return x.with_columns(
            narwhals.max_horizontal(
                (
                    narwhals.col(GWASLAB_95U_ODD_RATIO_COL).log()
                    - narwhals.col(GWASLAB_95L_ODD_RATIO_COL)
                )
                / (2 * 1.96),
                narwhals.lit(self.min_se),
            ).alias(GWASLAB_SE_COL)
        )
