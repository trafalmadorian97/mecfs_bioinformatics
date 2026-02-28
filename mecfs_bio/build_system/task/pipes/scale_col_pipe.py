import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class ScaleColPipe(DataProcessingPipe):
    col: str
    scale_factor: float

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.with_columns(
            (self.scale_factor * narwhals.col(self.col)).alias(self.col)
        )
