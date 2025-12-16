import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class UNNestPipe(DataProcessingPipe):
    col_to_unnest: str

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.explode(self.col_to_unnest)
