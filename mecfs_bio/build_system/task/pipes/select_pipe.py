import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen(slots=True)
class SelectColPipe(DataProcessingPipe):
    cols_to_select: list[str]

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        y = x.select(*self.cols_to_select)
        return y
