import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class HeadPipe(DataProcessingPipe):
    num_rows: int

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.head(self.num_rows)
