import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class CastPipe(DataProcessingPipe):
    target_column: str
    type: narwhals.dtypes.DType
    new_col_name: str

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.with_columns(
            narwhals.col(self.target_column).cast(self.type).alias(self.new_col_name)
        )
