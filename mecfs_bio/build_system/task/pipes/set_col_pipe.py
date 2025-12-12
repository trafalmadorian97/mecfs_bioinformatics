import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class SetColToConstantPipe(DataProcessingPipe):
    col_name: str
    constant: str | int | float

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.with_columns(narwhals.lit(self.constant).alias(self.col_name))
