import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class MaxWithConstantPipe(DataProcessingPipe):
    col: str
    value: float
    new_col_name: str

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.with_columns(
            narwhals.max_horizontal(
                narwhals.col(self.col), narwhals.lit(self.value)
            ).alias(self.new_col_name)
        )
