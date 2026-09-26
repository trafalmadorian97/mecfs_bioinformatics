import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen(slots=True)
class DropColPipe(DataProcessingPipe):
    cols_to_drop: list[str]

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        y = x.drop(*self.cols_to_drop)
        return y
